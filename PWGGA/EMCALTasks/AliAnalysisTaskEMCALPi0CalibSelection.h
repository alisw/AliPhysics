#ifndef ALIANALYSISTASKEMCALPI0CALIBSELECTION_H
#define ALIANALYSISTASKEMCALPI0CALIBSELECTION_H

// $Id$

// Root includes
class TH1F;
#include "TH2I.h"
#include "TObjArray.h"

// AliRoot includes
#include "AliAnalysisTaskSE.h"
class AliEMCALGeometry;
#include "AliEMCALGeoParams.h"
class AliEMCALRecoUtils;

class AliAnalysisTaskEMCALPi0CalibSelection : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskEMCALPi0CalibSelection(const char* name);
  virtual ~AliAnalysisTaskEMCALPi0CalibSelection();
    
  void    CorrectClusters();
  
  void    FillHistograms();
  
  void    InitGeometryMatrices();
  
  void    InitTemperatureCorrections();
  
  void    UserCreateOutputObjects();
  
  void    UserExec(Option_t * opt);
    
  void    PrintInfo();

  void    Terminate(Option_t* opt);
  
  void    GetMaxEnergyCellPosAndClusterPos(AliVCaloCells* cells, AliVCluster* clu, Int_t& iSM, Int_t& ieta, Int_t& iphi);
  
  // Analysis parameter setting
  
  void    SetPairDTimeCut(Float_t t)                     { fDTimeCut    = t          ; }
  void    SetClusterMinTime(Float_t tmin)                { fTimeMin     = tmin       ; }
  void    SetClusterMaxTime(Float_t tmax)                { fTimeMax     = tmax       ; }

  void    SetAsymmetryCut(Float_t asy)                   { fAsyCut      = asy        ; }
  
  void    SetClusterMinEnergy(Float_t emin)              { fEmin        = emin       ; }
  void    SetClusterMaxEnergy(Float_t emax)              { fEmax        = emax       ; }
  
  void    SetClusterLambda0Cuts(Float_t min, Float_t max){ fL0max       = max        ;
                                                           fL0min       = min        ; }
  void    SetClusterMinNCells(Int_t n)                   { fMinNCells   = n          ; }
  void    SetNCellsGroup(Int_t n)                        { fGroupNCells = n          ; }
  void    SetLogWeight(Float_t w)                        { fLogWeight   = w          ; }
  	  
  void    SetPairMinMassCut(Float_t min)                 { fInvMassCutMin = min      ; }
  void    SetPairMaxMassCut(Float_t max)                 { fInvMassCutMax = max      ; }
  
  void    SwitchOnSameSM()                               { fSameSM = kTRUE           ; }
  void    SwitchOffSameSM()                              { fSameSM = kFALSE          ; }
  
  void    UseFilteredEventAsInput()                      { fFilteredInput = kTRUE    ; }
  void    UseNormalEventAsInput()                        { fFilteredInput = kFALSE   ; }
  
  void    SetTriggerName(TString name)                   { fTriggerName = name       ; }

  //Geometry setters
  
  void    SetGeometryName(TString name)                  { fEMCALGeoName = name      ; }
  TString GeometryName() const                           { return fEMCALGeoName      ; }
  void    SwitchOnLoadOwnGeometryMatrices()              { fLoadMatrices = kTRUE     ; }
  void    SwitchOffLoadOwnGeometryMatrices()             { fLoadMatrices = kFALSE    ; }
  void    SetGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fMatrix[i]    = m         ; }

  void    SetOADBFilePath(TString path)                  { fOADBFilePath      = path ; }
  
  // Cluster recalculation
  
  void    SwitchOnClusterCorrection()                    { fCorrectClusters = kTRUE  ; }
  void    SwitchOffClusterCorrection()                   { fCorrectClusters = kFALSE ; }
  void    SetEMCALRecoUtils(AliEMCALRecoUtils * ru)      { fRecoUtils = ru           ; }
  AliEMCALRecoUtils* GetEMCALRecoUtils()    const        { return fRecoUtils         ; }
  
  void    SetInvariantMassHistoBinRange(Int_t nBins, Float_t minbin, Float_t maxbin){
                                        fNbins     = nBins ; fMinBin     = minbin ; fMaxBin     = maxbin ; }

  void    SetTimeHistoBinRange         (Int_t nBins, Float_t minbin, Float_t maxbin){
                                        fNTimeBins = nBins ; fMinTimeBin = minbin ; fMaxTimeBin = maxbin ; }

  
  // Mask clusters
  void    SetNMaskCellColumns(Int_t n) {
    if(n > fNMaskCellColumns){ delete [] fMaskCellColumns ; fMaskCellColumns = new Int_t[n] ; }
    fNMaskCellColumns = n ; }
  
  void    SetMaskCellColumn(Int_t ipos, Int_t icol) { if(ipos < fNMaskCellColumns) fMaskCellColumns[ipos] = icol            ;
                                                      else printf("Not set, position larger than allocated set size first") ; }
  
  Bool_t  MaskFrameCluster(Int_t iSM, Int_t ieta) const;
  
private:

  AliEMCALGeometry  * fEMCALGeo;         //! EMCAL geometry
  TGeoHMatrix       * fMatrix[AliEMCALGeoParams::fgkEMCALModules]; // Geometry matrices with alignments
  Bool_t              fLoadMatrices;    //  Matrices set from configuration, not get from geometry.root or from ESDs/AODs
  
  TString             fEMCALGeoName;    //  Name of geometry to use.
  TString             fTriggerName;     //  Trigger name must contain this name
    
  AliEMCALRecoUtils * fRecoUtils;       //  Access to reconstruction utilities
  TString             fOADBFilePath ;   //  Default path $ALICE_ROOT/OADB/EMCAL, if needed change
  Bool_t              fCorrectClusters; //  Correct clusters energy, position etc.


  TRefArray         * fCaloClustersArr; //! list of clusters
  AliVCaloCells     * fEMCALCells;      //! list of cells
  
  TList             * fCuts ;           //! List with analysis cuts
  TList             * fOutputContainer; //! histogram container
  Double_t            fVertex[3];       //! primary vertex
  Bool_t              fFilteredInput;   //  Read input produced with filter.

  // analysis cuts
  
  Float_t   fEmin;               // min. cluster energy (GeV)
  Float_t   fEmax;               // max. cluster energy (GeV)
  Float_t   fL0min;              // min. cluster L0
  Float_t   fL0max;              // max. cluster L0

  Float_t   fDTimeCut;           // Maximum difference between time of cluster pairs (ns)
  Float_t   fTimeMax;            // Maximum cluster time (ns)
  Float_t   fTimeMin;            // Minimum cluster time (ns)

  Float_t   fAsyCut;             // Asymmetry cut
  Int_t     fMinNCells;          // min. ncells in cluster
  Int_t     fGroupNCells;        // group n cells
  Float_t   fLogWeight;          // log weight used in cluster recalibration
  Bool_t    fSameSM;             // Combine clusters in channels on same SM

  Int_t     fNMaskCellColumns;   // Number of masked columns
  Int_t*    fMaskCellColumns;    //[fNMaskCellColumns] list of masked cell collumn
    
  Float_t   fInvMassCutMin;      // Min mass cut for clusters to fill time or other histograms
  Float_t   fInvMassCutMax;      // Mas mass cut for clusters to fill time or other histograms
  
  //Output histograms	and settings

  Int_t     fNbins;              // N       mass bins of invariant mass histograms
  Float_t   fMinBin;             // Minimum mass bins of invariant mass histograms
  Float_t   fMaxBin;             // Maximum mass bins of invariant mass histograms
  
  Int_t     fNTimeBins;          // N       time bins of invariant mass histograms
  Float_t   fMinTimeBin;         // Minimum time bins of invariant mass histograms
  Float_t   fMaxTimeBin;         // Maximum time bins of invariant mass histograms
  
  TH1F*     fHmpi0[AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows];//! two-cluster inv. mass assigned to each cell.

  TH2F*     fHmgg;                                                                 //! two-cluster inv.mass vs pt of pair
  TH2F*     fHmggDifferentSM;                                                      //! two-cluster inv.mass vs pt of pair, each cluster in different SM
  TH2F*     fHmggSM[AliEMCALGeoParams::fgkEMCALModules];                           //! two-cluster inv.mass per SM
  TH2F*     fHmggPairSameSectorSM[AliEMCALGeoParams::fgkEMCALModules/2];           //! two-cluster inv.mass per Pair
  TH2F*     fHmggPairSameSideSM  [AliEMCALGeoParams::fgkEMCALModules-2];           //! two-cluster inv.mass per Pair
  
  TH2F*     fHmggMaskFrame;                                                        //! two-cluster inv.mass vs pt of pair, mask clusters facing frames
  TH2F*     fHmggDifferentSMMaskFrame;                                             //! two-cluster inv.mass vs pt of pair, each cluster in different SM,mask clusters facing frames
  TH2F*     fHmggSMMaskFrame[AliEMCALGeoParams::fgkEMCALModules];                  //! two-cluster inv.mass per SM, mask clusters facing frames
  TH2F*     fHmggPairSameSectorSMMaskFrame[AliEMCALGeoParams::fgkEMCALModules/2];  //! two-cluster inv.mass per Pair, mask clusters facing frames
  TH2F*     fHmggPairSameSideSMMaskFrame  [AliEMCALGeoParams::fgkEMCALModules-2];  //! two-cluster inv.mass per Pair, mask clusters facing frames

  TH2F*     fHOpeningAngle;                                                        //! two-cluster opening angle vs pt of pair, with mass close to pi0
  TH2F*     fHOpeningAngleDifferentSM;                                             //! two-cluster opening angle vs pt of pair, each cluster in different SM, with mass close to pi0
  TH2F*     fHOpeningAngleSM[AliEMCALGeoParams::fgkEMCALModules];                  //! two-cluster opening angle vs pt per SM,with mass close to pi0
  TH2F*     fHOpeningAnglePairSM[AliEMCALGeoParams::fgkEMCALModules];              //! two-cluster opening angle vs pt per Pair,with mass close to pi0

  TH2F*     fHAsymmetry;                                                           //! two-cluster asymmetry vs pt of pair, with mass close to pi0
  TH2F*     fHAsymmetryDifferentSM;                                                //! two-cluster asymmetry vs pt of pair, each cluster in different SM, with mass close to pi0
  TH2F*     fHAsymmetrySM[AliEMCALGeoParams::fgkEMCALModules];                     //! two-cluster asymmetry vs pt per SM,with mass close to pi0
  TH2F*     fHAsymmetryPairSM[AliEMCALGeoParams::fgkEMCALModules];                 //! two-cluster asymmetry vs pt per Pair,with mass close to pi0
  
  TH2F*     fhTowerDecayPhotonHit[AliEMCALGeoParams::fgkEMCALModules] ;            //! Cells ordered in column/row for different module, number of times a decay photon hits
  TH2F*     fhTowerDecayPhotonEnergy[AliEMCALGeoParams::fgkEMCALModules] ;         //! Cells ordered in column/row for different module, accumulated energy in the tower by decay photons.
  TH2F*     fhTowerDecayPhotonAsymmetry[AliEMCALGeoParams::fgkEMCALModules] ;      //! Cells ordered in column/row for different module, accumulated asymmetry in the tower by decay photons.
  TH2F*     fhTowerDecayPhotonHitMaskFrame[AliEMCALGeoParams::fgkEMCALModules] ;   //! Cells ordered in column/row for different module, number of times a decay photon hits

  TH1I*     fhNEvents;                                                             //! Number of events counter histogram
 
  //Time
  TH2F*     fHTpi0[4];                                                             //! Time of cell under pi0 mass, for 4 bunch crossings
  TH2F*     fhClusterTime ;                                                        //! Timing of clusters vs energy
  TH2F*     fhClusterTimeSM[AliEMCALGeoParams::fgkEMCALModules] ;                  //! Timing of clusters vs energy per SM
  TH2F*     fhClusterPairDiffTime;                                                 //! Diference in time of clusters
  TH2F*     fhClusterPairDiffTimeSameSM[AliEMCALGeoParams::fgkEMCALModules];       //! Diference in time of clusters same SM
  TH2F*     fhClusterPairDiffTimeSameSector[AliEMCALGeoParams::fgkEMCALModules/2]; //! Diference in time of clusters same sector
  TH2F*     fhClusterPairDiffTimeSameSide[AliEMCALGeoParams::fgkEMCALModules-2];   //! Diference in time of clusters same side

  AliAnalysisTaskEMCALPi0CalibSelection(           const AliAnalysisTaskEMCALPi0CalibSelection&); 
  AliAnalysisTaskEMCALPi0CalibSelection& operator=(const AliAnalysisTaskEMCALPi0CalibSelection&); 
  
  ClassDef(AliAnalysisTaskEMCALPi0CalibSelection,20);

};

#endif //ALIANALYSISTASKEMCALPI0CALIBSELECTION_H
