#ifndef ALIANACALORIMETERQA_H
#define ALIANACALORIMETERQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

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

// --- Analysis system --- 
class AliVCaloCluster;
class AliVTrack;

#include "AliAnaPartCorrBaseClass.h"
 
class AliAnaCalorimeterQA : public AliAnaPartCorrBaseClass {
  
public: 
  AliAnaCalorimeterQA() ; // default ctor	
  virtual ~AliAnaCalorimeterQA() {;} //virtual dtor
private:
  AliAnaCalorimeterQA & operator = (const AliAnaCalorimeterQA & g) ;//cpy assignment
  AliAnaCalorimeterQA(const AliAnaCalorimeterQA & g) ; // cpy ctor
  
public:
  
  // General methods
  
  TObjString * GetAnalysisCuts();

  TList      * GetCreateOutputObjects();
  
  void         Init();
  
  void         InitParameters();
    
  void         MakeAnalysisFillHistograms() ;
  
  void         Print(const Option_t * opt) const;
    
  // Main methods
  
  void         ClusterHistograms(const TLorentzVector mom, Float_t *pos, 
                                 const Int_t nCaloCellsPerCluster, const Int_t nModule,
                                 const Bool_t matched, const AliVTrack* track, 
                                 const Int_t * labels, const Int_t nLabels);
 
  void         Correlate();

  void         FillCellPositionHistograms(const Int_t nCaloCellsPerCluster, const UShort_t * indexList, 
                                          const Float_t pos[3], const Float_t clEnergy);
  
  void         MCHistograms(const TLorentzVector mom, const Int_t pdg);
  
  // Setters and Getters

  
  Float_t      GetEMCALCellAmpMin()      const  { return fEMCALCellAmpMin    ; }
  void         SetEMCALCellAmpMin(Float_t amp)  { fEMCALCellAmpMin = amp     ; }
  
  Float_t      GetPHOSCellAmpMin()       const  { return fPHOSCellAmpMin     ; }
  void         SetPHOSCellAmpMin (Float_t amp)  { fPHOSCellAmpMin  = amp     ; }
    
  TString      GetCalorimeter()          const  { return fCalorimeter        ; }
  void         SetCalorimeter(TString calo)     { fCalorimeter = calo        ; }
    
  void         SetNumberOfModules(Int_t nmod)   { fNModules   = nmod         ; }
  
  Double_t     GetTimeCutMin()           const  { return fTimeCutMin         ; }
  Double_t     GetTimeCutMax()           const  { return fTimeCutMax         ; }
  void         SetTimeCut(Double_t min, Double_t max) {
                          fTimeCutMin = min ; fTimeCutMax = max              ; }
  
  // Histogram Switchs
  
  void SwitchOnFillAllPositionHistogram()       { fFillAllPosHisto  = kTRUE  ; }
  void SwitchOffFillAllPositionHistogram()      { fFillAllPosHisto  = kFALSE ; }
  
  void SwitchOnFillAllPositionHistogram2()      { fFillAllPosHisto2 = kTRUE  ; }
  void SwitchOffFillAllPositionHistogram2()     { fFillAllPosHisto2 = kFALSE ; }
  
  void SwitchOnFillAllTH12Histogram()           { fFillAllTH12      = kTRUE  ; }
  void SwitchOffFillAllTH12Histogram()          { fFillAllTH12      = kFALSE ; }
  
  void SwitchOnFillAllTH3Histogram()            { fFillAllTH3       = kTRUE  ; }
  void SwitchOffFillAllTH3Histogram()           { fFillAllTH3       = kFALSE ; }
  
  void SwitchOnFillAllTrackMatchingHistogram()  { fFillAllTMHisto   = kTRUE  ; }
  void SwitchOffFillAllTrackMatchingHistogram() { fFillAllTMHisto   = kFALSE ; }
  
  void SwitchOnFillAllPi0Histogram()            { fFillAllPi0Histo  = kTRUE  ; }
  void SwitchOffFillAllPi0Histogram()           { fFillAllPi0Histo  = kFALSE ; }

  void SwitchOnCorrelation()                    { fCorrelate        = kTRUE  ; }
  void SwitchOffCorrelation()                   { fCorrelate        = kFALSE ; }
    
 private:
  
  TString  fCalorimeter ;                     // Calorimeter selection
  Bool_t   fFillAllPosHisto;                  // Fill all the position related histograms 
  Bool_t   fFillAllPosHisto2;                 // Fill all the position related histograms 2
  Bool_t   fFillAllTH12 ;                     // Fill simple histograms which information is already in TH3 histograms
  Bool_t   fFillAllTH3 ;                      // Fill TH3 histograms
  Bool_t   fFillAllTMHisto ;                  // Fill track matching histograms
  Bool_t   fFillAllPi0Histo ;                 // Fill track matching histograms
  Bool_t   fCorrelate   ;                     // Correlate PHOS/EMCAL cells/clusters, also with V0 and track multiplicity
  Int_t    fNModules    ;                     // Number of EMCAL/PHOS modules
  Int_t    fNRCU        ;                     // Number of EMCAL/PHOS RCU 
  Int_t    fNMaxCols    ;                     // Number of EMCAL/PHOS rows 
  Int_t    fNMaxRows    ;                     // Number of EMCAL/PHOS columns
  Double_t fTimeCutMin  ;                     // Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;                     // Remove clusters/cells with time larger than this value, in ns
  Float_t  fEMCALCellAmpMin;                  // amplitude Threshold on emcal cells
  Float_t  fPHOSCellAmpMin ;                  // amplitude Threshold on phos cells
  
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
  TH2F *   fhNCellsPerClusterMIP;             //! N cells per cluster vs cluster energy vs eta of cluster, finer fixed pT bin for MIP search.
  TH2F *   fhNCellsPerClusterMIPCharged;      //! N cells per cluster vs cluster energy vs eta of cluster, finer fixed pT bin for MIP search, cluster matched with track.	
  
  TH2F *   fhNCellsvsClusterMaxCellDiffE0;    //! N cells per cluster vs cluster energy minus max cell, E < 2 GeV	
  TH2F *   fhNCellsvsClusterMaxCellDiffE2;    //! N cells per cluster vs cluster energy minus max cell, 2< E < 6	GeV
  TH2F *   fhNCellsvsClusterMaxCellDiffE6;    //! N cells per cluster vs cluster energy minus max cell, E > 6 GeV	
  
  TH1F *   fhNClusters;                       //! Number of clusters

  TH2F *   fhClusterTimeEnergy;               //! Cluster Time vs Energy 
  TH2F *   fhCellTimeSpreadRespectToCellMax;  //! Difference of the time of cell with maximum dep energy and the rest of cells
  TH2F *   fhClusterMaxCellDiffAverageTime;   //! Difference between cluster average time and time of cell with more energy
  TH2F *   fhClusterMaxCellDiffWeightTime;    //! Difference between cluster weighted average time and time of cell with more energy
  TH2F *   fhClusterDiffWeightAverTime;       //! Difference between cluster weighted average time and average time divided by weighted time
  TH2F *   fhClusterMaxCellDiffAverageNoMaxTime;//! Difference between cluster average time without max cell and time of cell with more energy
  TH2F *   fhClusterMaxCellDiffWeightNoMaxTime; //! Difference between cluster weighted average time without max cell and time of cell with more energy
  TH2F *   fhClusterNoMaxCellWeight;          //! energy weight of cells in cluster except maximum cell
  TH1F *   fhCellIdCellLargeTimeSpread;       //! Cells with large time respect to max (diff > 100 ns)
  TH2F *   fhClusterPairDiffTimeE;            //! Pair of clusters time difference vs E

  TH2F *   fhClusterMaxCellCloseCellRatio;    //! Ratio between max cell energy and cell energy of the same cluster 
  TH2F *   fhClusterMaxCellCloseCellDiff;     //! Difference between max cell energy and cell energy of the same cluster 

  TH2F *   fhClusterMaxCellDiff;              //! Difference between cluster energy and energy of cell with more energy, good clusters only
  TH2F *   fhClusterMaxCellDiffNoCut;         //! Difference between cluster energy and energy of cell with more energy, no bad cluster rejection
  
  TH2F *   fhLambda0vsClusterMaxCellDiffE0;   //! Lambda0 of bad cluster vs Fraction of energy of max cell for E < 2, no cut on bad clusters
  TH2F *   fhLambda0vsClusterMaxCellDiffE2;   //! Lambda0 of bad cluster vs Fraction of energy of max cell for E > 2, E < 6, no cut on bad clusters
  TH2F *   fhLambda0vsClusterMaxCellDiffE6;   //! Lambda0 of bad cluster vs Fraction of energy of max cell for E > 6, no cut on bad clusters
  
  TH1F *   fhBadClusterEnergy;                //! energy of bad cluster
  TH2F *   fhBadClusterTimeEnergy;            //! Time Max cell of bad cluster
  TH2F *   fhBadClusterPairDiffTimeE;         //! Pair of clusters time difference vs E, bad cluster
  TH2F *   fhBadClusterMaxCellCloseCellRatio; //! Ratio between max cell energy and cell energy of the same cluster for bad clusters 
  TH2F *   fhBadClusterMaxCellCloseCellDiff ; //! Difference between max cell energy and cell energy of the same cluster for bad clusters 
  TH2F *   fhBadClusterMaxCellDiff;           //! Difference between cluster energy and energy of cell with more energy
  TH2F *   fhBadClusterMaxCellDiffAverageTime;//! Difference between cluster average time and time of cell with more energy
  TH2F *   fhBadClusterMaxCellDiffWeightTime; //! Difference between cluster weighted average time and time of cell with more energy
  TH2F *   fhBadClusterDiffWeightAverTime;    //! Difference between cluster weighted average time and average time without max cell
  TH2F *   fhBadClusterMaxCellDiffAverageNoMaxTime;//! Difference between cluster average time without max cell and time of cell with more energy
  TH2F *   fhBadClusterMaxCellDiffWeightNoMaxTime; //! Difference between cluster weighted average time without max cell and time of cell with more energy
  TH2F *   fhBadClusterNoMaxCellWeight;       //! energy weight of cells in cluster except maximum cell
  TH2F *   fhBadCellTimeSpreadRespectToCellMax; //! Difference of the time of cell with maximum dep energy and the rest of cells for bad clusters

  TH2F *   fhBadClusterL0;                    //! Lambda0 for bad clusters
  TH2F *   fhBadClusterL1;                    //! Lambda1 for bad clusters
  TH2F *   fhBadClusterD;                     //! Dispersion for bad clusters

  // Cluster cell size
  TH2F *   fhDeltaIEtaDeltaIPhiE0[2];         // Difference between max cell index and farthest cell, eta vs phi, E < 2 GeV, with and without matching; 
  TH2F *   fhDeltaIEtaDeltaIPhiE2[2];         // Difference between max cell index and farthest cell, eta vs phi, 2 < E < 6 GeV, with and without matching; 
  TH2F *   fhDeltaIEtaDeltaIPhiE6[2];         // Difference between max cell index and farthest cell, eta vs phi, E > 6 GeV, with and without matching; 
  TH2F *   fhDeltaIA[2];                      // Cluster "asymmetry" in cell terms vs E, with and without matching
  TH2F *   fhDeltaIAL0[2];                    // Cluster "asymmetry" in cell units vs Lambda0    for E > 0.5 GeV, n cells in cluster > 3, with and without matching
  TH2F *   fhDeltaIAL1[2];                    // Cluster "asymmetry" in cell units vs Lambda1    for E > 0.5 GeV, n cells in cluster > 3, with and without matching
  TH2F *   fhDeltaIANCells[2] ;               // Cluster "asymmetry" in cell units vs number of cells in cluster for E > 0.5, with and without matching
  TH2F *   fhDeltaIAMC[4];                    // Cluster "asymmetry" in cell terms vs E, from MC photon, electron, conversion or hadron

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
  TH1F *   fhAmplitude;                       //! Amplitude measured in towers/crystals
  TH2F *   fhAmpId;                           //! Amplitude measured in towers/crystals vs id of tower.	
  TH3F *   fhEtaPhiAmp;                       //! eta vs phi vs amplitude, cells
  
  TH1F *   fhTime;                            //! Time measured in towers/crystals
  TH2F *   fhTimeId;                          //! Time vs Absolute cell Id
  TH2F *   fhTimeAmp;                         //! Time vs Amplitude 
  
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
  
  //Module histograms
  TH2F *   fhEMod  ;                          //! E distribution for different module, Reco
  TH2F *   fhNClustersMod ;                   //! Number of clusters for different module, Reco
  TH2F **  fhNCellsPerClusterMod ;            //! N cells per clusters different module, Reco
  TH2F **  fhNCellsPerClusterModNoCut ;       //! N cells per clusters different module, Reco, No cut
  TH2F *   fhNCellsMod ;                      //! Number of towers/crystals with signal different module, Reco
  TH2F *   fhGridCellsMod ;                   //! Cells ordered in column/row for different module, Reco
  TH2F *   fhGridCellsEMod ;                  //! Cells ordered in column/row for different module, weighted with energy, Reco
  TH2F *   fhGridCellsTimeMod ;               //! Cells ordered in column/row for different module, weighted with time, Reco
  TH2F **  fhTimeAmpPerRCU;                   //! Time vs Amplitude measured in towers/crystals different RCU
  TH2F **  fhIMMod;                           //! cluster pairs invariant mass, different module,
	
  //Pure MC

  enum mcTypes {mcPhoton = 0, mcPi0 = 1, mcEta = 2, mcElectron = 3, mcNeHadron = 4, mcChHadron = 5 };
  
  TH2F *   fhRecoMCE[6][2]  ;                 //! E   generated particle vs reconstructed E
  TH2F *   fhRecoMCPhi[6][2];                 //! phi generated particle vs reconstructed phi
  TH2F *   fhRecoMCEta[6][2];                 //! eta generated particle vs reconstructed Eta
  TH2F *   fhRecoMCDeltaE[6][2]  ;            //! Gen-Reco E    generated particle vs reconstructed E
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
//  TH2F *   fh1pOverE;                         //! p/E for track-cluster matches
//  TH1F *   fh1dR;                             //! distance between projected track and cluster
//  TH2F *   fh2EledEdx;                        //! dE/dx vs. momentum for electron candidates
//  TH2F *   fh2MatchdEdx;                      //! dE/dx vs. momentum for all matches
//	
//  TH2F *   fhMCEle1pOverE;                    //! p/E for track-cluster matches, MC electrons
//  TH1F *   fhMCEle1dR;                        //! distance between projected track and cluster, MC electrons
//  TH2F *   fhMCEle2MatchdEdx;                 //! dE/dx vs. momentum for all matches, MC electrons	
//	
//  TH2F *   fhMCChHad1pOverE;                  //! p/E for track-cluster matches, MC charged hadrons
//  TH1F *   fhMCChHad1dR;                      //! distance between projected track and cluster, MC charged hadrons
//  TH2F *   fhMCChHad2MatchdEdx;               //! dE/dx vs. momentum for all matches, MC charged
//	
//  TH2F *   fhMCNeutral1pOverE;                //! p/E for track-cluster matches, MC neutral
//  TH1F *   fhMCNeutral1dR;                    //! distance between projected track and cluster, MC neutral
//  TH2F *   fhMCNeutral2MatchdEdx;             //! dE/dx vs. momentum for all matches, MC neutral	
//	
//  TH2F *   fh1pOverER02;                      //! p/E for track-cluster matches, dR > 0.2	
//  TH2F *   fhMCEle1pOverER02;                 //! p/E for track-cluster matches, dR > 0.2, MC electrons
//  TH2F *   fhMCChHad1pOverER02;               //! p/E for track-cluster matches, dR > 0.2, MC charged hadrons
//  TH2F *   fhMCNeutral1pOverER02;             //! p/E for track-cluster matches, dR > 0.2, MC neutral
	
  ClassDef(AliAnaCalorimeterQA,19)
} ;


#endif //ALIANACALORIMETERQA_H



