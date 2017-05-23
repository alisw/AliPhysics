#ifndef ALIANACALORIMETERQA_H
#define ALIANACALORIMETERQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaCalorimeterQA
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Class for the Calorimeter QA analysis
///
/// Task filling basic histograms to check reconstructed data (real or MC) of calorimeters.
/// EMCal official QA executed with this task although it also works for PHOS.
/// The output of this task is also used for the tagging of channels as bad.
/// Output and conclusions of this task when executed on data can be
/// found in the following pages
///  * [Main CaloQA page](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CaloQA)
///  * [Web repository](http://aliqaemc.web.cern.ch/aliqaemc/)
///  * [QA trends](https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalQATrends)
///  * [Period by period](https://twiki.cern.ch/twiki/bin/view/ALICE/EMCALQAPeriodbyPeriod)
///  * [Run by run](https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalQARunByRun)
///  * [Bad channels](https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalQABadChannels)
///  * [Some documentation](https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalQAHowTo)
///
/// More information on the code can be found in this
/// [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this
/// [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaCalorimeterQA).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
///_________________________________________________________________________

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
    
  AliAnaCalorimeterQA() ;
    
  /// Virtual destructor. Not implemented.
  virtual ~AliAnaCalorimeterQA() { ; }
    
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
    
  void         ClusterAsymmetryHistograms(AliVCluster* clus, Int_t absIdMax, Bool_t goodCluster );
  
  void         ClusterHistograms(AliVCluster* cluster, const TObjArray *caloClusters,  AliVCaloCells * cells, 
                                 Int_t absIdMax, Double_t maxCellFraction, Float_t eCrossFrac, Double_t tmax);
  
  void         ClusterLoopHistograms(const TObjArray * clusters, AliVCaloCells * cells);
  
  Bool_t       ClusterMCHistograms(Bool_t matched, const Int_t * labels, Int_t nLabels, Int_t & pdg );

  void         ClusterMatchedWithTrackHistograms(AliVCluster* clus, Bool_t mcOK, Int_t pdg);

  void         Correlate();
        
  void         InvariantMassHistograms(Int_t iclus, Int_t nModule, const TObjArray* caloClusters, AliVCaloCells * cells);

  Bool_t       IsGoodCluster(Int_t absIdMax, Float_t m02, Int_t nCellsPerCluster, AliVCaloCells *cells);

  void         MCHistograms();  
  
  // Setters and getters
  
  Float_t      GetEMCALCellAmpMin()      const  { return fEMCALCellAmpMin    ; }
  void         SetEMCALCellAmpMin(Float_t amp)  { fEMCALCellAmpMin = amp     ; }
  
  Float_t      GetPHOSCellAmpMin()       const  { return fPHOSCellAmpMin     ; }
  void         SetPHOSCellAmpMin (Float_t amp)  { fPHOSCellAmpMin  = amp     ; }

  Float_t      GetEMCALM02Min()          const  { return fEMCALClusterM02Min ; }
  void         SetEMCALM02Min(Float_t m02)      { fEMCALClusterM02Min = m02  ; }
  
  Int_t        GetEMCALNCellsPerClusterMin() const  { return fEMCALClusterNCellMin ; }
  void         SetEMCALNCellsPerClusterMin(Int_t n) { fEMCALClusterNCellMin = n    ; }
  
  Int_t        GetPHOSNCellsPerClusterMin()  const  { return fPHOSClusterNCellMin  ; }
  void         SetPHOSNCellsPerClusterMin (Int_t n) { fPHOSClusterNCellMin  = n    ; }
  
  Float_t      GetInvMassMinECut()       const  { return fInvMassMinECut     ; }
  void         SetInvMassMinECut(Float_t cut)   { fInvMassMinECut = cut      ; }

  Float_t      GetInvMassMaxECut()       const  { return fInvMassMaxECut     ; }
  void         SetInvMassMaxECut(Float_t cut)   { fInvMassMaxECut = cut      ; }

  Float_t      GetInvMassMinM02Cut()     const  { return fInvMassMinM02Cut   ; }
  void         SetInvMassMinM02Cut(Float_t cut) { fInvMassMinM02Cut = cut    ; }
  
  Float_t      GetInvMassMaxM02Cut()     const  { return fInvMassMaxM02Cut   ; }
  void         SetInvMassMaxM02Cut(Float_t cut) { fInvMassMaxM02Cut = cut    ; }
  
  Float_t      GetInvMassMaxOpenAngle()  const  { return fInvMassMaxOpenAngle; }
  void         SetInvMassMaxOpenAngle(Float_t c){ fInvMassMaxOpenAngle = c   ; }
  
  Float_t      GetInvMassMaxTimeDifference()  const   { return fInvMassMaxTimeDifference; }
  void         SetInvMassMaxTimeDifference(Float_t c) { fInvMassMaxTimeDifference = c   ; }  
  
  Double_t     GetTimeCutMin()           const  { return fTimeCutMin         ; }
  Double_t     GetTimeCutMax()           const  { return fTimeCutMax         ; }
  void         SetTimeCut(Double_t min, Double_t max) {
                          fTimeCutMin = min ; fTimeCutMax = max              ; }
    
  void         SetNEBinCuts(Int_t nb)           { fNEBinCuts = nb            ; }
  void         SetEBinCutsAt(Int_t i, Float_t va) { if(i < 15) fEBinCuts[i] = va ; }

  
  // Histogram switchs
  
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

  void SwitchOnFillAllClusterHistogram()        { fFillAllClusterHistograms = kTRUE  ; }
  void SwitchOffFillAllClusterHistogram()       { fFillAllClusterHistograms = kFALSE ; }
  
  void SwitchOnFillAllCellHistogram()           { fFillAllCellHistograms  = kTRUE  ; }
  void SwitchOffFillAllCellHistogram()          { fFillAllCellHistograms  = kFALSE ; }
  
  void SwitchOnFillAllCellAbsIdHistogram()      { fFillAllCellAbsIdHistograms = kTRUE  ; }
  void SwitchOffFillAllCellAbsIdHistogram()     { fFillAllCellAbsIdHistograms = kFALSE ; }
  
  void SwitchOnFillInvMassOpAngleHistogram()    { fFillInvMassOpenAngle = kTRUE  ; }
  void SwitchOffFillInvMassOpAngleHistogram()   { fFillInvMassOpenAngle = kFALSE ; }  

  void SwitchOnFillPi0PairDiffTimeHistogram()   { fFillPi0PairDiffTime = kTRUE  ; }
  void SwitchOffFillPi0PairDiffTimeHistogram()  { fFillPi0PairDiffTime = kFALSE ; }

  void SwitchOnFillInvMassReducedEMCALHistogram()  { fFillInvMassInEMCALWithPHOSDCalAcc = kTRUE  ; }
  void SwitchOffFillInvMassReducedEMCALHistogram() { fFillInvMassInEMCALWithPHOSDCalAcc = kFALSE ; }

  void SwitchOnCorrelation()                    { fCorrelate        = kTRUE  ; }
  void SwitchOffCorrelation()                   { fCorrelate        = kFALSE ; }

  void SwitchOnStudyBadClusters()               { fStudyBadClusters = kTRUE  ; }
  void SwitchOffStudyBadClusters()              { fStudyBadClusters = kFALSE ; }

  void SwitchOnFillClusterCellMaxHisto()        { fFillClusterMaxCellHisto = kTRUE  ; }
  void SwitchOffFillClusterCellMaxHisto()       { fFillClusterMaxCellHisto = kFALSE ; }
  
  void SwitchOnAcceptanceHistoPerEBin()         { fFillEBinAcceptanceHisto = kTRUE  ; }
  void SwitchOffAcceptanceHistoPerEBin()        { fFillEBinAcceptanceHisto = kFALSE ; }
          
  void SetConstantTimeShift(Float_t shift)      { fConstantTimeShift     = shift  ; }
  
 private:
  
  // Switches
    
  Bool_t   fFillAllCellTimeHisto;               ///<  Fill all cell time histo
  Bool_t   fFillAllPosHisto;                    ///<  Fill all the position related histograms
  Bool_t   fFillAllPosHisto2;                   ///<  Fill all the position related histograms 2
  Bool_t   fFillAllTH3 ;                        ///<  Fill TH3 histograms
  Bool_t   fFillAllTMHisto ;                    ///<  Fill track matching histograms

  Bool_t   fFillClusterMaxCellHisto ;           ///<  Fill cluster cell max histograms
  Bool_t   fFillAllPi0Histo ;                   ///<  Fill invariant mass histograms
  Bool_t   fFillInvMassOpenAngle;               ///<  Fill opening angle histograms of cluster pairs, only if  fFillAllPi0Histo=kTRUE
  Bool_t   fFillPi0PairDiffTime;                ///<  Fill time difference histograms of cluster pairs in pi0 mass window, only if  fFillAllPi0Histo=kTRUE
  Bool_t   fFillInvMassInEMCALWithPHOSDCalAcc;  ///<  Fill invariant mass histograms of EMCal clusters in DCal and PHOS eta acceptance each, only if  fFillAllPi0Histo=kTRUE
  
  Bool_t   fFillEBinAcceptanceHisto;            ///<  Fill histograms with cluster eta-phi distribution and column-row cell, for different energy bins
  
  Bool_t   fFillAllClusterHistograms;           ///<  Fill all cluster related histograms
  Bool_t   fFillAllCellHistograms;              ///<  Fill all cell related histograms
  Bool_t   fFillAllCellAbsIdHistograms;         ///<  Fill all cell related histograms where one axis is the cell absId
  
  Bool_t   fCorrelate   ;                       ///<  Correlate PHOS/EMCAL cells/clusters, also with V0 and track multiplicity
  Bool_t   fStudyBadClusters;                   ///<  Study bad clusters not passing selection criteria (exotic, shower shape, n cells). 
  
  // Cuts
  
  Double_t fTimeCutMin  ;                       ///<  Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;                       ///<  Remove clusters/cells with time larger than this value, in ns
  Float_t  fCellAmpMin;                         ///<  Amplitude Threshold on calorimeter cells, set at execution time
  Float_t  fEMCALCellAmpMin;                    ///<  Amplitude Threshold on EMCal cells
  Float_t  fPHOSCellAmpMin ;                    ///<  Amplitude Threshold on PHOS cells
  Float_t  fEMCALClusterM02Min;                 ///<  Minimum M02 on EMCal clusters
  Int_t    fEMCALClusterNCellMin;               ///<  Minimum number of cells on EMCal clusters
  Int_t    fPHOSClusterNCellMin ;               ///<  Minimum number of cells on PHOS clusters
  Float_t  fEBinCuts[15] ;                      ///<  Energy bins cut 
  Int_t    fNEBinCuts;                          ///<  Number of energy bin cuts

  // Invariant mass analysis
  
  Float_t  fInvMassMinECut;                     ///<  Minimum energy cut value for clusters entering the invariant mass calculation
  Float_t  fInvMassMaxECut;                     ///<  Maximum energy cut value for clusters entering the invariant mass calculation
  Float_t  fInvMassMinM02Cut;                   ///<  Minimum M02 shower shape cut value for clusters entering the invariant mass calculation
  Float_t  fInvMassMaxM02Cut;                   ///<  Maximum M02 shower shape cut value for clusters entering the invariant mass calculation
  Float_t  fInvMassMaxOpenAngle;                ///<  Combine clusters within with a maximum opening angle between them. In radians.
  Float_t  fInvMassMaxTimeDifference;           ///<  Maximum difference between the time of the 2 clusters to be considered in invariant mass. In ns.
  
  TLorentzVector fClusterMomentum;              //!<! Cluster momentum, temporary container
  TLorentzVector fClusterMomentum2;             //!<! Cluster momentum, temporary container
  TLorentzVector fPrimaryMomentum;              //!<! Primary MC momentum, temporary container
  
  Float_t  fConstantTimeShift;                  ///<  Apply a 600 ns time shift in case of simulation, shift in ns.

  
  // Calorimeter Clusters
    
  TH1F *   fhE  ;                               //!<! E distribution, Reco
  TH1F *   fhPt ;                               //!<! pT distribution, Reco
  TH1F *   fhPhi;                               //!<! phi distribution, Reco
  TH1F *   fhEta;                               //!<! eta distribution, Reco
  TH2F *   fhEtaPhi;                            //!<! eta-phi distribution, Reco
  TH3F *   fhEtaPhiE  ;                         //!<! eta vs phi vs E, Reco
  TH1F *   fhECharged  ;                        //!<! E distribution, Reco, matched with track
  TH1F *   fhPtCharged ;                        //!<! pT distribution, Reco, matched with track
  TH1F *   fhPhiCharged;                        //!<! phi distribution, Reco, matched with track
  TH1F *   fhEtaCharged;                        //!<! eta-phi distribution, Reco, matched with track
  TH2F *   fhEtaPhiCharged;                     //!<! eta distribution, Reco, matched with track
  TH3F *   fhEtaPhiECharged;                    //!<! eta vs phi vs E, Reco, matched with track
    
  TH2F *   fhIM;                                //!<! Cluster pairs invariant mass vs pair pT, for EMCAL or PHOS pairs
  TH2F *   fhIMSame;                            //!<! Cluster pairs invariant mass vs pair pT, for EMCAL or PHOS pairs
  TH2F *   fhIMDiff;                            //!<! Cluster pairs invariant mass vs pair pT, for EMCAL or PHOS pairs
  TH2F *   fhIMDCAL;                            //!<! Cluster pairs invariant mass vs pair pT, for DCal pairs
  TH2F *   fhIMDCALSame;                        //!<! Cluster pairs invariant mass vs pair pT, for DCal pairs
  TH2F *   fhIMDCALDiff;                        //!<! Cluster pairs invariant mass vs pair pT, for DCal pairs
  TH2F *   fhIMDCALPHOS;                        //!<! Cluster pairs invariant mass vs pair pT, for DCal-PHOS pairs
  TH2F *   fhIMDCALPHOSSame;                    //!<! Cluster pairs invariant mass vs pair pT, for DCal-PHOS pairs
  TH2F *   fhIMEMCALPHOS;                       //!<! Cluster pairs invariant mass vs pair pT, for EMCAL(DCal eta acceptance)-EMCAL (PHOS eta acceptance) pairs
  TH2F *   fhIMEMCALPHOSSame;                   //!<! Cluster pairs invariant mass vs pair pT, for EMCAL(DCal eta acceptance)-EMCAL (PHOS eta acceptance) pairs
  
  TH2F *   fhAsym;                              //!<! Cluster pairs invariant mass vs pair pT
  
  TH2F*    fhOpAngle;                           //!<! Cluster pairs opening angle vs pair pT
  TH2F*    fhIMvsOpAngle;                       //!<! Cluster pairs opening angle vs mass  
  
  TH2F *   fhClusterPairDiffTimeEPi0Mass;       //!<! EMCal/PHOS Cluster time TOF difference, for pairs in 0.1 < mass < 0.18
  TH2F *   fhClusterPairDiffTimeEPi0MassSame;   //!<! EMCal/PHOS Cluster time TOF difference, for pairs in 0.1 < mass < 0.18, pairs in same Module
  
  TH2F *   fhClusterPairDiffTimeEPi0MassDCal;     //!<! DCal Cluster time TOF difference, for pairs in 0.1 < mass < 0.18
  TH2F *   fhClusterPairDiffTimeEPi0MassDCalSame; //!<! DCal Cluster time TOF difference, for pairs in 0.1 < mass < 0.18, pairs in same Module
  
  TH2F *   fhNCellsPerCluster;                  //!<! N cells per cluster vs cluster energy 
  TH2F *   fhNCellsPerClusterNoCut;             //!<! N cells per cluster vs cluster energy, before cuts 
 
  TH1F *   fhNClusters;                         //!<! Number of clusters

  TH2F *   fhClusterTimeEnergy;                 //!<! Cluster Time vs Energy
  TH2F *   fhCellTimeSpreadRespectToCellMax;    //!<! Difference of the time of cell with maximum dep energy and the rest of cells
  TH1F *   fhCellIdCellLargeTimeSpread;         //!<! Cells with large time respect to max (diff > 100 ns)
  TH2F *   fhClusterPairDiffTimeE;              //!<! Pair of clusters time difference vs E
  TH2F *   fhClusterPairDiffTimeESameMod;       //!<! Pair of clusters time difference vs E, in same Mod
  
  TH2F *   fhClusterMaxCellCloseCellRatio;      //!<! Ratio between max cell energy and cell energy of the same cluster
  TH2F *   fhClusterMaxCellCloseCellDiff;       //!<! Difference between max cell energy and cell energy of the same cluster
  TH2F *   fhClusterMaxCellDiff;                //!<! Difference between cluster energy and energy of cell with more energy, good clusters only
  TH2F *   fhClusterMaxCellDiffNoCut;           //!<! Difference between cluster energy and energy of cell with more energy, no bad cluster rejection

//TH2F *   fhClusterMaxCellDiffAverageTime;     //!<! Difference between cluster average time and time of cell with more energy
//TH2F *   fhClusterMaxCellDiffWeightedTime;    //!<! Difference between cluster weighted time and time of cell with more energy
  TH2F *   fhClusterMaxCellECross;              //!<! 1 - Energy in cross around max energy cell / max energy cell vs cluster energy, good clusters
  
//TH2F *   fhDispersion;                        //!<! Cluster Dispersion vs Energy
  TH2F *   fhLambda0;                           //!<! Cluster Lambda0 vs Energy
  TH2F *   fhLambda1;                           //!<! Cluster Lambda1 vs Energy  
  TH2F *   fhNLocMax;                           //!<! Cluster Number of local Maxima
  
  // Bad clusters histograms
  
  TH1F *   fhBadClusterEnergy;                  //!<! Energy of bad cluster
  TH2F *   fhBadClusterTimeEnergy;              //!<! Time Max cell of bad cluster
  TH2F *   fhBadClusterEtaPhi;                  //!<! Time Max cell of bad cluster
  TH2F *   fhBadClusterPairDiffTimeE;           //!<! Pair of clusters time difference vs E, bad cluster
  TH2F *   fhBadCellTimeSpreadRespectToCellMax; //!<! Difference of the time of cell with maximum dep energy and the rest of cells for bad clusters
  
  TH2F *   fhBadClusterMaxCellCloseCellRatio;   //!<! Ratio between max cell energy and cell energy of the same cluster for bad clusters
  TH2F *   fhBadClusterMaxCellCloseCellDiff ;   //!<! Difference between max cell energy and cell energy of the same cluster for bad clusters
  TH2F *   fhBadClusterMaxCellDiff;             //!<! Difference between cluster energy and energy of cell with more energy
  
//TH2F *   fhBadClusterMaxCellDiffAverageTime;  //!<! Difference between cluster average time and time of cell with more energy
//TH2F *   fhBadClusterMaxCellDiffWeightedTime; //!<! Difference between cluster weighted time and time of cell with more energy
  TH2F *   fhBadClusterMaxCellECross;           //!<! 1 - Energy in cross around max energy cell / max energy cell vs cluster energy, bad clusters

  TH2F *   fhBadClusterLambda0;                 //!<! Cluster Lambda0 vs Energy, clusters declared bad
  TH2F *   fhBadClusterLambda1;                 //!<! Cluster Lambda1 vs Energy, clusters declared bad

  // Cluster/cell Position
    
  TH2F *   fhRNCells ;                          //!<! R=sqrt(x^2+y^2) (cm) cluster distribution vs N cells in cluster
  TH2F *   fhXNCells ;                          //!<! X (cm) cluster distribution vs N cells in cluster
  TH2F *   fhYNCells ;                          //!<! Y (cm) cluster distribution vs N cells in cluster
  TH2F *   fhZNCells ;                          //!<! Z (cm) cluster distribution vs N cells in cluster
	
  TH2F *   fhRE ;                               //!<! R=sqrt(x^2+y^2) (cm) cluster distribution vs cluster energy
  TH2F *   fhXE ;                               //!<! X (cm) cluster distribution vs cluster energy
  TH2F *   fhYE ;                               //!<! Y (cm) cluster distribution vs cluster energy
  TH2F *   fhZE ;                               //!<! Z (cm) cluster distribution vs cluster energy
  TH3F *   fhXYZ;                               //!<! cluster X vs Y vs Z (cm)
	
  TH2F *   fhRCellE ;                           //!<! R=sqrt(x^2+y^2) (cm) cell distribution vs cell energy
  TH2F *   fhXCellE ;                           //!<! X (cm) cell distribution vs cell energy
  TH2F *   fhYCellE ;                           //!<! Y (cm) cell distribution vs cell energy
  TH2F *   fhZCellE ;                           //!<! Z (cm) cell distribution vs cell energy
  TH3F *   fhXYZCell;                           //!<! cell X vs Y vs Z (cm)
  
  TH2F *   fhDeltaCellClusterRNCells ;          //!<! R cluster - R cell distribution (cm) vs N cells in cluster
  TH2F *   fhDeltaCellClusterXNCells ;          //!<! X cluster - X cell distribution (cm) vs N cells in cluster
  TH2F *   fhDeltaCellClusterYNCells ;          //!<! Y cluster - Y cell distribution (cm) vs N cells in cluster
  TH2F *   fhDeltaCellClusterZNCells ;          //!<! Z cluster - Z cell distribution (cm) vs N cells in cluster
	
  TH2F *   fhDeltaCellClusterRE ;               //!<! R cluster - R cell distribution (cm) vs cluster energy
  TH2F *   fhDeltaCellClusterXE ;               //!<! X cluster - X cell distribution (cm) vs cluster energy
  TH2F *   fhDeltaCellClusterYE ;               //!<! Y cluster - Y cell distribution (cm) vs cluster energy
  TH2F *   fhDeltaCellClusterZE ;               //!<! Z cluster - Z cell distribution (cm) vs cluster energy
	
  // Calorimeter cells
    
  TH1F *   fhNCells;                            //!<! Number of towers/crystals with signal
  TH1F *   fhNCellsCutAmpMin;                   //!<! Number of towers/crystals with signal, with min amplitude
  TH1F *   fhAmplitude;                         //!<! Amplitude measured in towers/crystals
  TH2F *   fhAmpId;                             //!<! Amplitude measured in towers/crystals vs id of tower.
  TH3F *   fhEtaPhiAmpCell;                     //!<! eta vs phi vs amplitude, cells
  TH2F *   fhEtaPhiCell;                        //!<! eta vs phi, cells
   
  TH1F *   fhTime;                              //!<! Time measured in towers/crystals
//TH2F *   fhTimeVz;                            //!<! Time measured in towers/crystals vs vertex z component, for E > 0.5
  TH2F *   fhTimeId;                            //!<! Time vs Absolute cell Id
  TH2F *   fhTimeL1UnCorrId;                    //!<! Time (not corrected for L1 phase) vs Absolute cell Id
  TH2F *   fhTimeAmp;                           //!<! Time vs Amplitude
  TH2F *   fhTimePerSMPerBC[4];                 //!<! Time vs SM number for BC%4=0,1,2,3
  
  TH2F *   fhAmpIdLowGain;                      //!<! Amplitude measured in towers/crystals vs id of tower, low gain towers
  TH2F *   fhTimeIdLowGain;                     //!<! Time vs Absolute cell Id, low gain
  TH2F *   fhTimeAmpLowGain;                    //!<! Time vs Amplitude, low gain

  TH2F *   fhCellECross;                        //!<! 1 - Energy in cross around cell /  cell energy
  
  // Calorimeters Correlation
  
  TH2F *   fhEMCALPHOSCorrNClusters;            //!<! EMCAL vs PHOS, number of clusters
  TH2F *   fhEMCALPHOSCorrEClusters;            //!<! EMCAL vs PHOS, total measured cluster energy
  TH2F *   fhEMCALPHOSCorrNCells;               //!<! EMCAL vs PHOS, number of cells
  TH2F *   fhEMCALPHOSCorrECells;               //!<! EMCAL vs PHOS, total measured cell energy

  TH2F *   fhEMCALDCALCorrNClusters;            //!<! EMCAL vs DCAL, number of clusters
  TH2F *   fhEMCALDCALCorrEClusters;            //!<! EMCAL vs DCAL, total measured cluster energy
  TH2F *   fhEMCALDCALCorrNCells;               //!<! EMCAL vs DCAL, number of cells
  TH2F *   fhEMCALDCALCorrECells;               //!<! EMCAL vs DCAL, total measured cell energy

  TH2F *   fhDCALPHOSCorrNClusters;             //!<! DCAL vs PHOS, number of clusters
  TH2F *   fhDCALPHOSCorrEClusters;             //!<! DCAL vs PHOS, total measured cluster energy
  TH2F *   fhDCALPHOSCorrNCells;                //!<! DCAL vs PHOS, number of cells
  TH2F *   fhDCALPHOSCorrECells;                //!<! DCAL vs PHOS, total measured cell energy
  
  // V0 Correlation
  
  TH2F *   fhCaloV0SCorrNClusters;              //!<! Calo vs V0 signal , number of clusters
  TH2F *   fhCaloV0SCorrEClusters;              //!<! Calo vs V0 signal, total measured cluster energy
  TH2F *   fhCaloV0SCorrNCells;                 //!<! Calo vs V0 signal, number of cells
  TH2F *   fhCaloV0SCorrECells;                 //!<! Calo vs V0 signal,  total measured cell energy
  TH2F *   fhCaloV0MCorrNClusters;              //!<! Calo vs V0 multiplicity , number of clusters
  TH2F *   fhCaloV0MCorrEClusters;              //!<! Calo vs V0 multiplicity, total measured cluster energy
  TH2F *   fhCaloV0MCorrNCells;                 //!<! Calo vs V0 multiplicity, number of cells
  TH2F *   fhCaloV0MCorrECells;                 //!<! Calo vs V0 multiplicity,  total measured cell energy
  
  // Track Correlation
  
  TH2F *   fhCaloTrackMCorrNClusters;           //!<! Calo vs Track Multiplicity, number of clusters
  TH2F *   fhCaloTrackMCorrEClusters;           //!<! Calo vs Track Multiplicity, total measured cluster energy
  TH2F *   fhCaloTrackMCorrNCells;              //!<! Calo vs V0 Track Multiplicity, number of cells
  TH2F *   fhCaloTrackMCorrECells;              //!<! Calo vs V0 Track Multipliticy,  total measured cell energy
  
  // Centrality
  
  TH2F *   fhCaloCenNClusters;                  //!<! Calo vs centrality, number of clusters
  TH2F *   fhCaloCenEClusters;                  //!<! Calo vs centrality, total measured cluster energy
  TH2F *   fhCaloCenNCells;                     //!<! Calo vs centrality, number of cells
  TH2F *   fhCaloCenECells;                     //!<! Calo vs centrality,  total measured cell energy

  // Event plane
  
  TH2F *   fhCaloEvPNClusters;                  //!<! Calo vs event plane angle, number of clusters
  TH2F *   fhCaloEvPEClusters;                  //!<! Calo vs event plane angle, total measured cluster energy
  TH2F *   fhCaloEvPNCells;                     //!<! Calo vs event plane angle, number of cells
  TH2F *   fhCaloEvPECells;                     //!<! Calo vs event plane angle,  total measured cell energy
  
  // Module histograms
  
  TH2F *   fhEMod  ;                            //!<! Cluster E distribution for different module, Reco
  TH2F *   fhAmpMod ;                           //!<! Cell amplitude distribution for different module, Reco
  TH2F *   fhEWeirdMod  ;                       //!<! Cluster E distribution for different module, very large E, Reco
  TH2F *   fhAmpWeirdMod ;                      //!<! Cell amplitude distribution for different module, very large Amp, Reco
  TH2F *   fhTimeMod ;                          //!<! Cell time distribution for different module, Reco
  TH2F *   fhNClustersMod ;                     //!<! Number of clusters for different module, Reco
  TH2F *   fhNCellsMod ;                        //!<! Number of towers/crystals with signal for different module, Reco  
  TH2F *   fhSumClustersEnergyMod ;             //!<! Sum of clusters  energy for different module, Reco
  TH2F *   fhSumCellsAmpMod ;                   //!<! Sum of towers/crystals signal for different module, Reco
  TH2F **  fhNCellsSumAmpPerMod ;               //!<! N cells vs sum of amplitude in different modules, Reco
  TH2F **  fhNClustersSumEnergyPerMod ;         //!<! N clusters vs sum of energies in different module, Reco
  TH2F **  fhNCellsPerClusterMod ;              //!<! N cells per clusters different module, Reco
  TH2F **  fhNCellsPerClusterModNoCut ;         //!<! N cells per clusters different module, Reco, No cut
  TH2F *   fhNCellsPerClusterWeirdMod ;         //!<! N cells per clusters different module, Reco, ridiculously large energy
  TH2F *   fhNCellsPerClusterWeirdModNoCut ;    //!<! N cells per clusters different module, Reco, No cut, ridiculously large energy
  TH2F *   fhGridCells ;                        //!<! Cells ordered in column/row for different module, Reco
  TH2F *   fhGridCellsE ;                       //!<! Cells ordered in column/row for different module, weighted with energy, Reco
  TH2F *   fhGridCellsTime ;                    //!<! Cells ordered in column/row for different module, weighted with time, Reco
  TH2F *   fhGridCellsLowGain ;                 //!<! Cells ordered in column/row for different module, Reco, low gain
  TH2F *   fhGridCellsELowGain ;                //!<! Cells ordered in column/row for different module, weighted with energy, Reco, low gain
  TH2F *   fhGridCellsTimeLowGain ;             //!<! Cells ordered in column/row for different module, weighted with time, Reco, low gain
  TH2F **  fhTimeAmpPerRCU;                     //!<! Time vs Amplitude measured in towers/crystals different RCU
  TH2F **  fhIMMod;                             //!<! cluster pairs invariant mass, different module,
	
  // Pure MC histograms

  /// Enumerator with indeces for MC histograms array indicating the particle type generating the cluster
  enum mcTypes { kmcPhoton   = 0, kmcPi0        = 1, kmcEta = 2,
                 kmcElectron = 3, kmcPhotonConv = 4,
                 kmcNeHadron = 5, kmcChHadron   = 6             };
  
  TH2F *   fhRecoMCE[7][2]  ;                   //!<! E   generated particle vs reconstructed E
  TH2F *   fhRecoMCPhi[7][2] ;                  //!<! phi generated particle vs reconstructed phi
  TH2F *   fhRecoMCEta[7][2] ;                  //!<! eta generated particle vs reconstructed Eta
  TH2F *   fhRecoMCDeltaE[7][2]  ;              //!<! Gen-Reco E    generated particle vs reconstructed E
  TH2F *   fhRecoMCRatioE[7][2]  ;              //!<! Reco/Gen E    generated particle vs reconstructed E
  TH2F *   fhRecoMCDeltaPhi[7][2];              //!<! Gen-Reco phi  generated particle vs reconstructed E
  TH2F *   fhRecoMCDeltaEta[7][2];              //!<! Gen-Reco eta  generated particle vs reconstructed E
  
  TH1F *   fhGenMCE [4]     ;                   //!<! pt of primary particle
  TH1F *   fhGenMCPt[4]     ;                   //!<! pt of primary particle
  TH2F *   fhGenMCEtaPhi[4] ;                   //!<! eta vs phi of primary particle
  TH1F *   fhGenMCAccE [4]     ;                //!<! pt of primary particle, in acceptance
  TH1F *   fhGenMCAccPt[4]     ;                //!<! pt of primary particle, in acceptance
  TH2F *   fhGenMCAccEtaPhi[4] ;                //!<! eta vs phi of primary particle, in acceptance
  
  TH2F *   fhEMVxyz    ;                        //!<! Electromagnetic particle production vertex
  TH2F *   fhEMR       ;                        //!<! Electromagnetic distance to vertex vs rec energy
  TH2F *   fhHaVxyz    ;                        //!<! Hadron production vertex
  TH2F *   fhHaR       ;                        //!<! Hadron distance to vertex vs rec energy
	
  // Histograms for MC track-matching
  
  TH2F *   fh1EOverP;                           //!<! p/E for track-cluster matches
  TH2F *   fh2dR;                               //!<! distance between projected track and cluster (eta-phi units)
  TH2F *   fh2EledEdx;                          //!<! dE/dx vs. momentum for electron candidates
  TH2F *   fh2MatchdEdx;                        //!<! dE/dx vs. momentum for all matches
  TH2F *   fh1EOverPR02;                        //!<! p/E for track-cluster matches, dR < 0.2
  TH2F *   fh1EleEOverP;                        //!<! p/E for track-cluster matches, dR < 0.2, 60 < dEdx < 100

  TH2F *   fh1EOverPMod;                        //!<! p/E for track-cluster matches, per SM
  TH2F *   fh2dRMod;                            //!<! distance between projected track and cluster (eta-phi units), per SM
  TH2F *   fh2EledEdxMod;                       //!<! dE/dx for electron candidates, per SM
  TH2F *   fh2MatchdEdxMod;                     //!<! dE/dx for all matches, per SM
  TH2F *   fh1EOverPR02Mod;                     //!<! p/E for track-cluster matches, dR < 0.2, per SM
  TH2F *   fh1EleEOverPMod;                     //!<! p/E for track-cluster matches, dR < 0.2, 60 < dEdx < 100, per SM

  TH2F *   fhMCEle1EOverP;                      //!<! p/E for track-cluster matches, MC electrons
  TH1F *   fhMCEle1dR;                          //!<! distance between projected track and cluster, MC electrons
  TH2F *   fhMCEle2MatchdEdx;                   //!<! dE/dx vs. momentum for all matches, MC electrons
	
  TH2F *   fhMCChHad1EOverP;                    //!<! p/E for track-cluster matches, MC charged hadrons
  TH1F *   fhMCChHad1dR;                        //!<! distance between projected track and cluster, MC charged hadrons
  TH2F *   fhMCChHad2MatchdEdx;                 //!<! dE/dx vs. momentum for all matches, MC charged
	
  TH2F *   fhMCNeutral1EOverP;                  //!<! p/E for track-cluster matches, MC neutral
  TH1F *   fhMCNeutral1dR;                      //!<! Distance between projected track and cluster, MC neutral
  TH2F *   fhMCNeutral2MatchdEdx;               //!<! dE/dx vs. momentum for all matches, MC neutral
	
  TH2F *   fhMCEle1EOverPR02;                   //!<! p/E for track-cluster matches, dR < 0.2, MC electrons
  TH2F *   fhMCChHad1EOverPR02;                 //!<! p/E for track-cluster matches, dR < 0.2, MC charged hadrons
  TH2F *   fhMCNeutral1EOverPR02;               //!<! p/E for track-cluster matches, dR < 0.2, MC neutral

  TH2F *   fhMCEle1EleEOverP;                   //!<! p/E for track-cluster matches, dR < 0.2, 60 < dEdx < 100, MC electrons
  TH2F *   fhMCChHad1EleEOverP;                 //!<! p/E for track-cluster matches, dR < 0.2, 60 < dEdx < 100, MC charged hadrons
  TH2F *   fhMCNeutral1EleEOverP;               //!<! p/E for track-cluster matches, dR < 0.2, 60 < dEdx < 100, MC neutral

  TH2F *   fhTrackMatchedDEtaNeg;               //!<! Eta distance between track and cluster vs cluster E, after and before photon cuts
  TH2F *   fhTrackMatchedDPhiNeg;               //!<! Phi distance between track and cluster vs cluster E, after and before photon cuts
  TH2F *   fhTrackMatchedDEtaDPhiNeg;           //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV, after and before
  
  TH2F *   fhTrackMatchedDEtaPos;               //!<! Eta distance between track and cluster vs cluster E, after and before photon cuts
  TH2F *   fhTrackMatchedDPhiPos;               //!<! Phi distance between track and cluster vs cluster E, after and before photon cuts
  TH2F *   fhTrackMatchedDEtaDPhiPos;           //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV, after and before

  TH2F *   fhTrackMatchedDEtaNegMod;            //!<! Eta distance between negative track and cluster vs module for E > 0.5 GeV
  TH2F *   fhTrackMatchedDPhiNegMod;            //!<! Phi distance between negative track and cluster vs module for E > 0.5 GeV
  TH2F *   fhTrackMatchedDEtaPosMod;            //!<! Eta distance between positive track and cluster vs module for E > 0.5 GeV
  TH2F *   fhTrackMatchedDPhiPosMod;            //!<! Phi distance between positive track and cluster vs module for E > 0.5 GeV
  
  TH2F *   fhEBinClusterEtaPhi[14] ;            //!<! Eta-Phi location of cluster in different energy bins.
  TH2F *   fhEBinClusterColRow[14] ;            //!<! Column and row location of cluster max E cell in different energy bins.
  TH2F *   fhEBinCellColRow   [14] ;            //!<! Column and row location of cell in different energy bins.
  
  /// Copy constructor not implemented.
  AliAnaCalorimeterQA & operator = (const AliAnaCalorimeterQA & qa) ;
    
  /// Assignment operator not implemented.
  AliAnaCalorimeterQA(              const AliAnaCalorimeterQA & qa) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaCalorimeterQA,38) ;
  /// \endcond

} ;

#endif //ALIANACALORIMETERQA_H



