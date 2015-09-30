#ifndef ALIANALYSISTASKEMCALPI0CALIBSELECTION_H
#define ALIANALYSISTASKEMCALPI0CALIBSELECTION_H

//---------------------------------------------------------------------------
/// \class AliAnalysisTaskEMCALPi0CalibSelection
/// \brief This task provides the input for the EMCal energy calibration with pi0 invariant mass analysis per channel
///
/// Fill histograms (one per cell) with two-cluster invariant mass in EMCal
/// using calibration coefficients of the previous iteration.
/// Histogram for a given cell is filled if the most of the energy of the 2 clusters
/// is deposited in this cell and the other cluster could be anywhere in EMCAL.
/// Several cuts are possible to apply to improve the invariant mass distribution
/// per channel:
/// * Asymmetry cut
/// * Mininum and maximum energy of the clusters
/// * Shower shape
/// * Difference in time of the clusters
/// * Restrict cluster pairs to same super module
///
/// This analysis needs to be executed in several iterations, in each one the mean mass peak position
/// is obtained and applied as calibration factor for the next iteration.
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalCalibrationPage).
///
/// Adapted from PHOS tasks developped by Boris Polishchuk.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//---------------------------------------------------------------------------

// Root includes
class TH1F;
#include "TH2I.h"
#include "TObjArray.h"
#include "TLorentzVector.h"

// AliRoot includes
#include "AliAnalysisTaskSE.h"
class AliEMCALGeometry;
#include "AliEMCALGeoParams.h"
class AliEMCALRecoUtils;

class AliAnalysisTaskEMCALPi0CalibSelection : public AliAnalysisTaskSE
{
    
public:

  AliAnalysisTaskEMCALPi0CalibSelection();
  
  AliAnalysisTaskEMCALPi0CalibSelection(const char* name);
    
  virtual ~AliAnalysisTaskEMCALPi0CalibSelection();
    
  void    CorrectClusters();
  
  void    FillHistograms();
  
  void    InitEnergyCalibrationFactors();

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

  // Geometry setters
  
  void    SetGeometryName(TString name)                  { fEMCALGeoName = name      ; }
  
  TString GeometryName() const                           { return fEMCALGeoName      ; }
  
  void    SwitchOnLoadOwnGeometryMatrices()              { fLoadMatrices = kTRUE     ; }
  
  void    SwitchOffLoadOwnGeometryMatrices()             { fLoadMatrices = kFALSE    ; }
  
  void    SetGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fMatrix[i]    = m         ; }

  void    SetOADBFilePath(TString path)                  { fOADBFilePath = path      ; }
  
  void    SetCalibrationFilePath( TString path )         { fCalibFilePath = path     ; }
  
  // Cluster recalculation
  
  void    SwitchOnClusterCorrection()                    { fCorrectClusters = kTRUE  ; }
 
  void    SwitchOffClusterCorrection()                   { fCorrectClusters = kFALSE ; }

  void    SwitchOnRecalculatePosition()                  { fRecalPosition   = kTRUE  ; }
  
  void    SwitchOffRecalculatePosition()                 { fRecalPosition   = kFALSE ; }

  void    SetEMCALRecoUtils(AliEMCALRecoUtils * ru)      { fRecoUtils = ru           ; }
  
  AliEMCALRecoUtils* GetEMCALRecoUtils()    const        { return fRecoUtils         ; }
  
  void    SetInvariantMassHistoBinRange(Int_t nBins, Float_t minbin, Float_t maxbin){
                                        fNbins     = nBins ; fMinBin     = minbin ; fMaxBin     = maxbin ; }

  void    SetTimeHistoBinRange         (Int_t nBins, Float_t minbin, Float_t maxbin){
                                        fNTimeBins = nBins ; fMinTimeBin = minbin ; fMaxTimeBin = maxbin ; }

  
  void    SetImportGeometryFromFile(Bool_t import, TString path = ""){
                                                           fImportGeometryFromFile = import ;
                                                           fImportGeometryFilePath = path   ; } 

  // Mask clusters
  
  void    SetNMaskCellColumns(Int_t n) ;
  
  void    SetMaskCellColumn(Int_t ipos, Int_t icol) ;
  
  Bool_t  MaskFrameCluster(Int_t iSM, Int_t ieta) const;
  
private:

  AliEMCALGeometry  * fEMCALGeo;         //!<! EMCAL geometry pointer.
    
  ///< Geometry matrices with alignments.
  TGeoHMatrix       * fMatrix[AliEMCALGeoParams::fgkEMCALModules];
    
  Bool_t              fLoadMatrices;     ///<  Matrices set from configuration, not get from geometry.root or from ESDs/AODs.
    
  
  TString             fEMCALGeoName;     ///<  Name of geometry to use.
    
  TString             fTriggerName;      ///<  Trigger name must contain this name.
    
  AliEMCALRecoUtils * fRecoUtils;        ///<  Access to reconstruction utilities.
    
  TString             fOADBFilePath ;    ///<  Default path $ALICE_PHYSICS/OADB/EMCAL, if needed change.
  
  TString             fCalibFilePath;    ///< Full path with file with energy calibration factors per channel from previous iteration.
    
  Bool_t              fCorrectClusters;  ///<  Correct clusters energy, position etc.

  Bool_t              fRecalPosition;    ///<  Switch on/off cluster position calculation, in case alignment matrices are not available.
  
  TRefArray         * fCaloClustersArr;  //!<! List of clusters.
    
  AliVCaloCells     * fEMCALCells;       //!<! List of cells.
  
  TList             * fCuts ;            //!<! List with analysis cuts.
    
  TList             * fOutputContainer;  //!<! Histogram container.
    
  Double_t            fVertex[3];        //!<! Primary vertex.
    
  Bool_t              fFilteredInput;    ///<  Read input produced with filter.

  Bool_t              fImportGeometryFromFile; ///<  Import geometry settings in geometry.root file.
  
  TString             fImportGeometryFilePath; ///<  Path fo geometry.root file.
  
  // Analysis cuts
  
  Float_t             fEmin;             ///<  Minimum cluster energy (GeV).
  Float_t             fEmax;             ///<  Maximum cluster energy (GeV).
    
  Float_t             fL0min;            ///<  Minimum cluster L0.
  Float_t             fL0max;            ///<  Maximum cluster L0.

  Float_t             fDTimeCut;         ///<  Maximum difference between time of cluster pairs (ns).
  Float_t             fTimeMax;          ///<  Maximum cluster time (ns).
  Float_t             fTimeMin;          ///<  Minimum cluster time (ns).

  Float_t             fAsyCut;           ///<  Asymmetry cut.
    
  Int_t               fMinNCells;        ///<  Minimum ncells in cluster.
  Int_t               fGroupNCells;      ///<  Group n cells.
  
  Float_t             fLogWeight;        ///<  Logarithmic weight used in cluster recalibration.
    
  Bool_t              fSameSM;           ///<  Combine clusters in channels on same SM.

  Int_t               fNMaskCellColumns; ///<  Number of masked columns.
    
  ///< List the masked columns.
  Int_t*              fMaskCellColumns;  //[fNMaskCellColumns]
    
  Float_t             fInvMassCutMin;    ///<  Minimum mass cut for clusters to fill time or other histograms.
  Float_t             fInvMassCutMax;    ///< Maximum mass cut for clusters to fill time or other histograms.
  
  // Output histograms and settings

  Int_t               fNbins;            ///<  N       mass bins of invariant mass histograms.
  Float_t             fMinBin;           ///<  Minimum mass bins of invariant mass histograms.
  Float_t             fMaxBin;           ///<  Maximum mass bins of invariant mass histograms.
  
  Int_t               fNTimeBins;        ///<  N       time bins of invariant mass histograms.
  Float_t             fMinTimeBin;       ///<  Minimum time bins of invariant mass histograms.
  Float_t             fMaxTimeBin;       ///<  Maximum time bins of invariant mass histograms.
  
  // Temporal TLorentzVectors, avoir recreation per event 
    
  TLorentzVector      fMomentum1 ;       ///<  Cluster kinematics, temporal
  TLorentzVector      fMomentum2 ;       ///<  Cluster kinematics, temporal
  TLorentzVector      fMomentum12;       ///<  Cluster pair kinematics, temporal
    
  // Histograms

  ///< Two-cluster invariant mass assigned to each cell.
  TH1F*     fHmpi0[AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; //!

  TH2F*     fHmgg;                                                                 //!<! Two-cluster invariant mass vs pt of pair.
  TH2F*     fHmggDifferentSM;                                                      //!<! Two-cluster invariant mass vs pt of pair, each cluster in different SM.
  TH2F*     fHmggSM[AliEMCALGeoParams::fgkEMCALModules];                           //!<! Two-cluster invariant mass per SM.
  TH2F*     fHmggPairSameSectorSM[AliEMCALGeoParams::fgkEMCALModules/2];           //!<! Two-cluster invariant mass per Pair.
  TH2F*     fHmggPairSameSideSM  [AliEMCALGeoParams::fgkEMCALModules-2];           //!<! Two-cluster invariant mass per Pair.
  
  TH2F*     fHmggMaskFrame;                                                        //!<! Two-cluster invariant mass vs pt of pair, mask clusters facing frames.
  TH2F*     fHmggDifferentSMMaskFrame;                                             //!<! Two-cluster invariant mass vs pt of pair, each cluster in different SM,mask clusters facing frames.
  TH2F*     fHmggSMMaskFrame[AliEMCALGeoParams::fgkEMCALModules];                  //!<! Two-cluster invariant mass per SM, mask clusters facing frames.
  TH2F*     fHmggPairSameSectorSMMaskFrame[AliEMCALGeoParams::fgkEMCALModules/2];  //!<! Two-cluster invariant mass per Pair, mask clusters facing frames.
  TH2F*     fHmggPairSameSideSMMaskFrame  [AliEMCALGeoParams::fgkEMCALModules-2];  //!<! Two-cluster invariant mass per Pair, mask clusters facing frames.

  TH2F*     fHOpeningAngle;                                                        //!<! Two-cluster opening angle vs pt of pair, with mass close to pi0.
  TH2F*     fHOpeningAngleDifferentSM;                                             //!<! Two-cluster opening angle vs pt of pair, each cluster in different SM, with mass close to pi0.
  TH2F*     fHOpeningAngleSM[AliEMCALGeoParams::fgkEMCALModules];                  //!<! Two-cluster opening angle vs pt per SM,with mass close to pi0.
  TH2F*     fHOpeningAnglePairSM[AliEMCALGeoParams::fgkEMCALModules];              //!<! Two-cluster opening angle vs pt per Pair,with mass close to pi0.

  TH2F*     fHAsymmetry;                                                           //!<! Two-cluster asymmetry vs pt of pair, with mass close to pi0.
  TH2F*     fHAsymmetryDifferentSM;                                                //!<! Two-cluster asymmetry vs pt of pair, each cluster in different SM, with mass close to pi0.
  TH2F*     fHAsymmetrySM[AliEMCALGeoParams::fgkEMCALModules];                     //!<! Two-cluster asymmetry vs pt per SM,with mass close to pi0.
  TH2F*     fHAsymmetryPairSM[AliEMCALGeoParams::fgkEMCALModules];                 //!<! Two-cluster asymmetry vs pt per Pair,with mass close to pi0.
  
  TH2F*     fhTowerDecayPhotonHit[AliEMCALGeoParams::fgkEMCALModules] ;            //!<! Cells ordered in column/row for different module, number of times a decay photon hits.
  TH2F*     fhTowerDecayPhotonEnergy[AliEMCALGeoParams::fgkEMCALModules] ;         //!<! Cells ordered in column/row for different module, accumulated energy in the tower by decay photons.
  TH2F*     fhTowerDecayPhotonAsymmetry[AliEMCALGeoParams::fgkEMCALModules] ;      //!<! Cells ordered in column/row for different module, accumulated asymmetry in the tower by decay photons.
  TH2F*     fhTowerDecayPhotonHitMaskFrame[AliEMCALGeoParams::fgkEMCALModules] ;   //!<! Cells ordered in column/row for different module, number of times a decay photon hits.

  TH1I*     fhNEvents;                                                             //!<! Number of events counter histogram.
 
  // Cluster time histograms
  TH2F*     fHTpi0[4];                                                             //!<! Time of cell under pi0 mass, for 4 bunch crossings.
  TH2F*     fhClusterTime ;                                                        //!<! Timing of clusters vs energy.
  TH2F*     fhClusterTimeSM[AliEMCALGeoParams::fgkEMCALModules] ;                  //!<! Timing of clusters vs energy per SM.
  TH2F*     fhClusterPairDiffTime;                                                 //!<! Diference in time of clusters.
  TH2F*     fhClusterPairDiffTimeSameSM[AliEMCALGeoParams::fgkEMCALModules];       //!<! Diference in time of clusters same SM.
  TH2F*     fhClusterPairDiffTimeSameSector[AliEMCALGeoParams::fgkEMCALModules/2]; //!<! Diference in time of clusters same sector.
  TH2F*     fhClusterPairDiffTimeSameSide[AliEMCALGeoParams::fgkEMCALModules-2];   //!<! Diference in time of clusters same side.

  /// Copy constructor not implemented.
  AliAnalysisTaskEMCALPi0CalibSelection(           const AliAnalysisTaskEMCALPi0CalibSelection&) ;
    
  /// Assignment operator not implemented.
  AliAnalysisTaskEMCALPi0CalibSelection& operator=(const AliAnalysisTaskEMCALPi0CalibSelection&) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEMCALPi0CalibSelection,23) ;
  /// \endcond

};

#endif //ALIANALYSISTASKEMCALPI0CALIBSELECTION_H
