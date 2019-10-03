#ifndef ALIANALYSISTASKEMCALCLUSTERIZE_H
#define ALIANALYSISTASKEMCALCLUSTERIZE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliAnalysisTaskEMCALClusterize
/// \ingroup EMCALPerformance 
/// \brief Reclusterize EMCal clusters, put them in a new branch for other following analysis
///
/// This analysis provides a new list of clusters to be used in other analysis running right after this task. The clusters
/// are recalibrated, bad channels removed, track-matching recalculated. Clusters are put in a new branch.
/// Optionally, new clusters branch will be stored in an output AOD file with other additionnal information.
///
/// Adapted from analysis class from Deepa Thomas.
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// Root
class TTree;
class TClonesArray;
#include <TRandom3.h>

// EMCAL
class AliEMCALGeometry;
class AliEMCALClusterizer;
class AliEMCALAfterBurnerUF;
class AliEMCALRecPoint;
class AliAODCaloCluster;
class AliCentrality;
class AliMultSelection;
class AliVCaloCells;

#include "AliEMCALRecParam.h"
#include "AliEMCALRecoUtils.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALClusterize : public AliAnalysisTaskSE {
    
 public:
    
  AliAnalysisTaskEMCALClusterize();
  AliAnalysisTaskEMCALClusterize(const char *name);
  virtual ~AliAnalysisTaskEMCALClusterize();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Init();
  virtual void   LocalInit()                                   { Init()                        ; }
  void           PrintParam();
  
  // Event methods, settings
  
  Bool_t         AcceptCell(Int_t absID, Bool_t badmap = kTRUE);
  Bool_t         AcceptEventEMCAL();
  void           SwitchOnSelectEMCALEvent()                    { fSelectEMCALEvent   = kTRUE   ; }
  void           SwitchOffSelectEMCALEvent()                   { fSelectEMCALEvent   = kFALSE  ; }
  void           SetEMCALEnergyCut(Float_t cut)                { fEMCALEnergyCut     = cut     ; }
  void           SetEMCALNcellsCut(Int_t cut)                  { fEMCALNcellsCut     = cut     ; }

  void           SwitchOnInputAODFilter()                      { fInputFromFilter    = kTRUE   ; }
  void           SwitchOffInputAODFilter()                     { fInputFromFilter    = kFALSE  ; }
  
  void           CheckAndGetEvent();
  
  Bool_t         IsExoticEvent();
  void           SwitchOnExoticEventsRemoval()                  { fRemoveExoticEvents= kTRUE   ; }
  void           SwitchOffExoticEventsRemoval()                 { fRemoveExoticEvents= kFALSE  ; } 
  
  Bool_t         IsLEDEvent(const Int_t run);
  void           SwitchOnLEDEventsRemoval()                     { fRemoveLEDEvents   = kTRUE   ; }
  void           SwitchOffLEDEventsRemoval()                    { fRemoveLEDEvents   = kFALSE  ; } 
  
  // OCDB, avoid acessing OCDB!
    
  void           AccessOCDB();
  void           SwitchOnAccessOCDB()                           { fAccessOCDB       = kTRUE    ; }
  void           SwitchOffAccessOCDB()                          { fAccessOCDB       = kFALSE   ; } 
  void           SetOCDBPath(const char *path)                  { fOCDBpath         = path     ; }
  
  // Geometry methods
    
  void           InitGeometry();
  void           SetGeometryName(TString name)                  { fGeomName = name             ; }
  TString        GeometryName()                          const  { return fGeomName             ; }  
  void           SwitchOnLoadOwnGeometryMatrices()              { fLoadGeomMatrices = kTRUE    ; }
  void           SwitchOffLoadOwnGeometryMatrices()             { fLoadGeomMatrices = kFALSE   ; } 
  void           SetGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fGeomMatrix[i]    = m        ; }

  void           SetImportGeometryFromFile(Bool_t  im, 
                                           TString pa = "")     { fImportGeometryFromFile = im ; 
                                                                  fImportGeometryFilePath = pa ; }    
  // Outout AOD branch methods
    
  void           SetAODBranchName(TString name)                 { fOutputAODBranchName = name  ; }
  void           SetAODCellsName(TString name)                  { fOutputAODCellsName  = name  ; }
  void           SetInputCaloCellsName(TString name)            { fInputCaloCellsName  = name  ; }
  void           FillAODFile(Bool_t yesno)                      { fFillAODFile         = yesno ; }
  void           FillAODCaloCells();
  void           FillAODHeader();
  void           SwitchOnFillAODHeader()                        { fFillAODHeader     = kTRUE   ; }
  void           SwitchOffFillAODHeader()                       { fFillAODHeader     = kFALSE  ; } 
  void           SwitchOnFillAODCaloCells()                     { fFillAODCaloCells  = kTRUE   ; }
  void           SwitchOffFillAODCaloCells()                    { fFillAODCaloCells  = kFALSE  ; } 

  void           SwitchOnRecalibrateWithClusterTime()           { fRecalibrateWithClusterTime  = kTRUE   ; }
  void           SwitchOffRecalibrateWithClusterTime()          { fRecalibrateWithClusterTime  = kFALSE  ; }
  
  // Algorithms settings
  
  AliEMCALRecParam * GetRecParam()                              { if(!fRecParam)  fRecParam  = new AliEMCALRecParam  ;
                                                                  return fRecParam             ; }
  
  AliEMCALRecoUtils* GetRecoUtils()                             { if(!fRecoUtils) fRecoUtils = new AliEMCALRecoUtils ;  
                                                                  return fRecoUtils            ; }
  void          ConfigureEMCALRecoUtils(Bool_t  bMC    = kFALSE, Bool_t  bExotic= kTRUE, Bool_t  bNonLin= kFALSE,  
                                        Bool_t  bRecalE= kTRUE , Bool_t  bBad   = kTRUE, Bool_t  bRecalT= kTRUE, Int_t   debug  = -1);
  
  void           InitClusterization();
  void           ClusterizeCells();
  void           ClusterUnfolding();
  void           JustUnfold(Bool_t yesno)                       { fJustUnfold        = yesno   ; }
  void           UpdateCells();
  
  void           SetConfigFileName(TString name)                { fConfigName        = name    ; }
  void           SetMaxEvent(Int_t max)                         { fMaxEvent          = max     ; }
  void           SetMinEvent(Int_t max)                         { fMinEvent          = max     ; }
  
  void           SwitchOnTrackMatching()                        { fDoTrackMatching   = kTRUE   ; }
  void           SwitchOffTrackMatching()                       { fDoTrackMatching   = kFALSE  ; } 

  void           SwitchOnUpdateCell()                           { fUpdateCell        = kTRUE   ; } 
  void           SwitchOffUpdateCell()                          { fUpdateCell        = kFALSE  ; }  

  // Cell selection after unfolding
    
  void           SwitchOnCellEnergySelection()                  { fSelectCell        = kTRUE   ; }
  void           SwitchOffCellEnergySelection()                 { fSelectCell        = kFALSE  ; } 
  void           SetCellCuts(Float_t e, Float_t frac)           { fSelectCellMinE    = e       ; 
                                                                  fSelectCellMinFrac = frac    ; }
  void           SetRejectBelowThreshold(Bool_t reject)         { fRejectBelowThreshold =reject ; }
    
  // OADB options settings
  
  void           AccessOADB() ;
  
  TString        GetPass()    ;
  
  void           SwitchOnEMCALOADB()                            { fAccessOADB        = kTRUE   ; }
  void           SwitchOffEMCALOADB()                           { fAccessOADB        = kFALSE  ; }
    
  void           SetOADBFilePath(TString path)                  { fOADBFilePath      = path    ; }
  
  void           SetConstantTimeShift(Float_t shift)            { fConstantTimeShift = shift   ; }

  // Centrality selection
  
  AliCentrality* GetCentrality()                          const { return InputEvent()->GetCentrality() ; } 
  AliMultSelection* GetMultSelCen()                       const { return (AliMultSelection * ) InputEvent()->FindListObject("MultSelection") ; }
  void           SwitchOnAliCentrality ()                       { fUseAliCentrality  = kTRUE  ; }
  void           SwitchOffAliCentrality()                       { fUseAliCentrality  = kFALSE ; }
  
  void           SetCentralityClass(TString name)               { fCentralityClass   = name            ; }
  TString        GetCentralityClass()                     const { return fCentralityClass              ; }
  Float_t        GetEventCentrality()                     const ;
  void           SetCentralityBin(Int_t min, Int_t max) //Set the centrality bin to select the event. If used, then need to get percentile
                                                                { fCentralityBin[0]=min ; fCentralityBin[1]=max ; }
  Float_t        GetCentralityBin(Int_t i)                const { if(i < 0 || i > 1) return -1 ; 
                                                                  else               return fCentralityBin[i]   ; }
  
  // MC label properly assignation methods
  
  void           RemapMCLabelForAODs(Int_t &label);
  void           SwitchOnRemapMCLabelForAODs()                  { fRemapMCLabelForAODs       = kTRUE   ; }
  void           SwitchOffRemapMCLabelForAODs()                 { fRemapMCLabelForAODs       = kFALSE  ; }

  void           SetClustersMCLabelFrom2SelectedLabels(AliEMCALRecPoint* recPoint, AliAODCaloCluster *clus) ;
  void           SetClustersMCLabelFromOriginalClusters(AliAODCaloCluster * clus) ;
  
  void           SwitchOnUseClusterMCLabelForCell(Int_t opt = 0){ fSetCellMCLabelFromCluster = opt     ; }
  void           SwitchOffUseClusterMCLabelForCell()            { fSetCellMCLabelFromCluster = 0       ; }

  void           SwitchOnUseMCEdepFracLabelForCell()            { fSetCellMCLabelFromEdepFrac = kTRUE  ;  
                                                                   fSetCellMCLabelFromCluster = 0      ; }
  void           SwitchOffUseMCEdepFracLabelForCell()           { fSetCellMCLabelFromEdepFrac = kFALSE ; }
  
  //-----------------------------------------
  // T-Card correlation emulation, do on MC
  
  void           MakeCellTCardCorrelation() ;
  void           CalculateInducedEnergyInTCardCell(Int_t absId, Int_t absIdRef, Int_t sm, Float_t ampRef, Int_t cellCase) ;
  void           AddNewTCardInducedCellsToDigit() ;
  
  /// Activate T-Card cells correlation, 
  /// \param conservEnergy activate cluster energy conservation, not by default
  void           SwitchOnTCardCorrelation(Bool_t conservEnergy = kFALSE)  { fTCardCorrEmulation = kTRUE  ; fTCardCorrClusEnerConserv = conservEnergy ; }   
  
  /// De-activate T-Card cells correlation, 
  void           SwitchOffTCardCorrelation()                              { fTCardCorrEmulation = kFALSE ; fTCardCorrClusEnerConserv = kFALSE        ; }      

  
  /// Constant energy lost by max energy cell in one of T-Card cells, same for all SM
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossConstant(Float_t ud, Float_t udlr, Float_t lr, Float_t sec) { 
    for(Int_t ism = 0; ism < 22; ism++) {
      fTCardCorrInduceEner[0][ism] = ud; fTCardCorrInduceEner[1][ism] = udlr;  
      fTCardCorrInduceEner[2][ism] = lr; fTCardCorrInduceEner[3][ism] = sec ; } } 
  
  /// Fraction of energy lost by max energy cell in one of T-Card cells, same for all SM
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossFraction(Float_t ud, Float_t udlr, Float_t lr, Float_t sec) { 
    for(Int_t ism = 0; ism < 22; ism++) {
      fTCardCorrInduceEnerFrac[0][ism] = ud; fTCardCorrInduceEnerFrac[1][ism] = udlr;  
      fTCardCorrInduceEnerFrac[2][ism] = lr; fTCardCorrInduceEnerFrac[3][ism] = sec ; } } 

  /// Slope parameter of fraction of energy lost by max energy cell in one of T-Card cells, same for all SM
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossFractionP1(Float_t ud, Float_t udlr, Float_t lr, Float_t sec) { 
    for(Int_t ism = 0; ism < 22; ism++) {
      fTCardCorrInduceEnerFracP1[0][ism] = ud; fTCardCorrInduceEnerFracP1[1][ism] = udlr;  
      fTCardCorrInduceEnerFracP1[2][ism] = lr; fTCardCorrInduceEnerFracP1[3][ism] = sec ; } }

  /// Constant energy lost by max energy cell in one of T-Card cells, per SM
  /// \param sm super module index
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossConstantPerSM(Int_t sm, Float_t ud, Float_t udlr, Float_t lr, Float_t sec) { 
    if ( sm < 22 && sm >= 0 ) {
      fTCardCorrInduceEner[0][sm] = ud; fTCardCorrInduceEner[1][sm] = udlr;  
      fTCardCorrInduceEner[2][sm] = lr; fTCardCorrInduceEner[3][sm] = sec ; } } 
  
  /// Fraction of energy lost by max energy cell in one of T-Card cells, per SM
  /// \param sm super module index
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossFractionPerSM(Int_t sm, Float_t ud, Float_t udlr, Float_t lr, Float_t sec) { 
    if ( sm < 22 && sm >= 0 ) {
      fTCardCorrInduceEnerFrac[0][sm] = ud; fTCardCorrInduceEnerFrac[1][sm] = udlr;  
      fTCardCorrInduceEnerFrac[2][sm] = lr; fTCardCorrInduceEnerFrac[3][sm] = sec ; } } 
  
  /// Slope parameter of fraction of energy lost by max energy cell in one of T-Card cells, per SM
  /// \param sm super module index
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossFractionP1PerSM(Int_t sm, Float_t ud, Float_t udlr, Float_t lr, Float_t sec) { 
    if ( sm < 22 && sm >= 0 ) {
      fTCardCorrInduceEnerFracP1[0][sm] = ud; fTCardCorrInduceEnerFracP1[1][sm] = udlr;  
      fTCardCorrInduceEnerFracP1[2][sm] = lr; fTCardCorrInduceEnerFracP1[3][sm] = sec ; } }

  /// Fraction of energy lost by max energy cell in one of T-Card cells, width of random gaussian, same for all SM
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossFractionWidth(Float_t ud, Float_t udlr, Float_t lr, Float_t sec) {
    for(Int_t ism = 0; ism < 22; ism++) {
      fTCardCorrInduceEnerFracWidth[0][ism] = ud; fTCardCorrInduceEnerFracWidth[1][ism] = udlr;  
      fTCardCorrInduceEnerFracWidth[2][ism] = lr; fTCardCorrInduceEnerFracWidth[3][ism] = sec ; } }

  /// Fraction of energy lost by max energy cell in one of T-Card cells, width of random gaussian, per SM
  /// \param sm super module index
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossFractionWidthPerSM(Int_t sm, Float_t ud, Float_t udlr, Float_t lr, Float_t sec) {
    if ( sm < 22 && sm >= 0 ) {
      fTCardCorrInduceEnerFracWidth[0][sm] = ud; fTCardCorrInduceEnerFracWidth[1][sm] = udlr;  
      fTCardCorrInduceEnerFracWidth[2][sm] = lr; fTCardCorrInduceEnerFracWidth[3][sm] = sec ; } }

  /// Maximum induced energy fraction when linear dependency is set, per SM number
  /// \param max maximum fraction
  /// \param sm  super-module number
  void           SetInducedEnergyLossMaximumFractionPerSM(Float_t max, Int_t sm) { 
    if ( sm < 22 && sm >= 0 ) fTCardCorrInduceEnerFracMax[sm] = max ; }  
  
  /// Minimum induced energy fraction when linear dependency is set, per SM number
  /// \param min minimum fraction
  /// \param sm  super-module number
  void           SetInducedEnergyLossMinimumFractionPerSM(Float_t min, Int_t sm) { 
    if ( sm < 22 && sm >= 0 ) fTCardCorrInduceEnerFracMin[sm] = min ; }  
  
  /// Maximum induced energy fraction when linear dependency is set, same for all SM
  /// \param max maximum fraction
  void           SetInducedEnergyLossMaximumFraction(Float_t max) { 
    for(Int_t ism = 0; ism < 22; ism++) fTCardCorrInduceEnerFracMax[ism] = max ; }  
  
  /// Minimum induced energy fraction when linear dependency is set, same for all SM
  /// \param min minimum fraction
  void           SetInducedEnergyLossMinimumFraction(Float_t min) { 
    for(Int_t ism = 0; ism < 22; ism++) fTCardCorrInduceEnerFracMin[ism] = min ; }  
  
  /// fraction of times max cell energy correlates with cross cells, different for each super-module
  /// \param prob probability per event, from 0 to 1
  /// \param sm   probability assigned to this super-module number
  void           SetInducedEnergyLossProbabilityPerSM(Float_t prob, Int_t sm) { 
    if ( sm < 22 && sm >= 0 ) fTCardCorrInduceEnerProb[sm] = prob ; }  
  
  void           SwitchOnRandomizeTCardInducedEnergy()          { fRandomizeTCard = kTRUE   ; } 
  void           SwitchOffRandomizeTCardInducedEnergy()         { fRandomizeTCard = kFALSE  ; }  

  void           SetInducedTCardMinimumCellEnergy(Float_t mi)   { fTCardCorrMinAmp     = mi ; }
  void           SetInducedTCardMaximum(Float_t ma)             { fTCardCorrMaxInduced = ma ; }
  void           SetInducedTCardMinimum(Float_t mi)             { fTCardCorrMinInduced = mi ; }
  void           SetInducedTCardMaximumLowE(Float_t ma)         { fTCardCorrMaxInducedLowE = ma ; }
  
  void           PrintTCardParam();

  void     SwitchUseMergedBCs(Bool_t doUseMergedBC)     { fDoMergedBCs     = doUseMergedBC; }

  void     SetUse1DRecalibration(Bool_t use1D)     { fLoad1DRecalibFactors     = use1D; }
  
  //------------------------------------------
  
private:
    
  virtual void   FillCaloClusterInEvent();
  
  virtual void   RecPoints2Clusters();
  
  virtual void   ResetArrays();
    
  AliVEvent             *fEvent;                   ///<  Event 
  
  // Geometry
  AliEMCALGeometry      *fGeom;                    ///<  EMCAL geometry
  TString                fGeomName;                ///<  Name of geometry to use.
  TGeoHMatrix           *fGeomMatrix[22];          ///<  Geometry matrices with alignments
  Bool_t                 fGeomMatrixSet;           ///<  Set geometry matrices only once, for the first event.         
  Bool_t                 fLoadGeomMatrices;        ///<  Matrices set from configuration, not get from geometry.root or from ESDs/AODs

  // OCDB
  TString                fOCDBpath;                ///<  Path with OCDB location
  Bool_t                 fAccessOCDB;              ///<  Need to access info from OCDB (not really)   

  // Temporal arrays
  TClonesArray          *fDigitsArr;               //!<! Digits array
  TObjArray             *fClusterArr;              //!<! Recpoints array
  TObjArray             *fCaloClusterArr;          //!<! CaloClusters array
  AliVCaloCells         *fCaloCells;               //!<! CaloCells container

  // Clusterizers
  AliEMCALRecParam      *fRecParam;                ///<  Reconstruction parameters container
  AliEMCALClusterizer   *fClusterizer;             //!<! EMCAL clusterizer
  AliEMCALAfterBurnerUF *fUnfolder;                //!<! Unfolding procedure
  Bool_t                 fJustUnfold;              ///<  Just unfold, do not recluster
  
  // AOD
  TClonesArray          *fOutputAODBranch;         //!<! AOD Branch with output clusters
  TString                fOutputAODBranchName;     ///<  New of output clusters AOD branch
  AliAODCaloCells       *fOutputAODCells;          //!<! AOD Branch with output cells
  TString                fOutputAODCellsName;      ///<  New of output cells AOD branch name
  TString                fInputCaloCellsName;      ///<  Input cells branch name, if different from default branch
  Bool_t                 fOutputAODBranchSet ;     ///<  Set the AOD clusters branch in the input event once
  Bool_t                 fFillAODFile;             ///<  Fill the output AOD file with the new clusters, 
                                                   ///<  if not they will be only available for the event they were generated
  Bool_t                 fFillAODHeader;           ///<  Copy header to standard branch
  Bool_t                 fFillAODCaloCells;        ///<  Copy calocells to standard branch

  Int_t                  fRun;                     ///<  run number
  
  AliEMCALRecoUtils*     fRecoUtils;               ///<  Access to factorized reconstruction algorithms
  TString                fConfigName;              ///<  Name of analysis configuration file
    
  static const Int_t fgkNEMCalCells = 17664;       ///< Total number of cells in the calorimeter, 10*48*24 (EMCal) + 4*48*8 (EMCal/DCal 1/3) + 6*32*24 (DCal)
  
  Int_t                  fOrgClusterCellId[fgkNEMCalCells]; ///<  Array ID of cluster to wich the cell belongs in unmodified clusters
  Int_t                  fCellLabels      [fgkNEMCalCells]; ///<  Array with MC label to be passed to digit.
  Int_t                  fCellSecondLabels[fgkNEMCalCells]; ///<  Array with Second MC label to be passed to digit.
  Double_t               fCellTime        [fgkNEMCalCells]; ///<  Array with cluster time to be passed to digit in case of AODs
  Float_t                fCellMatchdEta   [fgkNEMCalCells]; ///<  Array with cluster-track dPhi
  Float_t                fCellMatchdPhi   [fgkNEMCalCells]; ///<  Array with cluster-track dEta

  Bool_t                 fRecalibrateWithClusterTime;       ///<  Use fCellTime to store time of cells in cluster
  
  Int_t                  fMaxEvent;                ///<  Set a maximum event number, for testing
  Int_t                  fMinEvent;                ///<  Set a minimum event number, for testing
  
  Bool_t                 fDoTrackMatching;         ///<  On/Off the matching recalculation to speed up analysis in PbPb
  Bool_t                 fUpdateCell;              ///<  On/Off the upate of the CaloCells container
  Bool_t                 fSelectCell;              ///<  Reject cells from cluster if energy is too low and recalculate position/energy and other
  Float_t                fSelectCellMinE;          ///<  Min energy cell threshold, after unfolding
  Float_t                fSelectCellMinFrac;       ///<  Min fraction of cell energy after unfolding cut
  Bool_t                 fRejectBelowThreshold;    ///<  split (false-default) or reject (true) cell energy below threshold after UF
  Bool_t                 fRemoveLEDEvents;         ///<  Remove LED events, use only for LHC11a 
  Bool_t                 fRemoveExoticEvents;      ///<  Remove exotic events
  
  Bool_t                 fImportGeometryFromFile;  ///<  Import geometry settings in geometry.root file
  TString                fImportGeometryFilePath;  ///<  path fo geometry.root file

  Bool_t                 fOADBSet ;                ///<  AODB parameters already set
  Bool_t                 fAccessOADB ;             ///<  Get calibration from OADB for EMCAL
  TString                fOADBFilePath ;           ///<  Default path $ALICE_PHYSICS/OADB/EMCAL, if needed change
  Float_t                fConstantTimeShift;       ///<  Apply a 600 ns time shift in case of simulation, shift in ns.

  // Centrality
  TString                fCentralityClass;         ///<  Name of selected centrality class     
  Float_t                fCentralityBin[2];        ///<  Minimum and maximum value of the centrality for the analysis
  Bool_t                 fUseAliCentrality;        ///<  Use the centrality estimator from AliCentrality or AliMultSelection

  //  Event selection with some signal in EMCAL
  Bool_t                 fSelectEMCALEvent;        ///<   Process the event if there is some high energy cluster.
  Float_t                fEMCALEnergyCut;          ///<   At least an EMCAL cluster with this energy in the event.
  Int_t                  fEMCALNcellsCut;          ///<   At least an EMCAL cluster with fNCellsCut cells over fEnergyCut.

  ///<  Use cluster MC label as cell label:
  ///<   * 0 - get the MC label stored in cells
  ///<   * 1 - select 2 most likely labels
  ///<   * 2 - get the original clusters, add all the MC labels 
  ///< Options 1 and 2 useful for any reclusterization with output V1 clusters and similar clusterization thresholds as original cluster,
  ///< if original is 50 MeV cell E cut and new is 100 MeV, this does not work well.
  Int_t                  fSetCellMCLabelFromCluster;
  
  ///< For MC generated with aliroot > v5-07-21, check the EDep information 
  ///< stored in ESDs/AODs to set the cell MC labels
  Bool_t                 fSetCellMCLabelFromEdepFrac;  
    
  Bool_t                 fRemapMCLabelForAODs ;    ///<  Remap AOD cells MC label. Needed in old AOD productions.

  Bool_t                 fInputFromFilter ;        ///<  Get the input from AODs from the filter.
    
  
  // T-Card correlation emulation, do on MC
  Bool_t                fTCardCorrEmulation;       ///< Activate T-Card cells energy correlation
  Bool_t                fTCardCorrClusEnerConserv; ///< When making correlation, subtract from the reference cell the induced energy on the neighbour cells
  Float_t               fTCardCorrCellsEner[fgkNEMCalCells]; ///<  Array with induced cell energy in T-Card neighbour cells
  Bool_t                fTCardCorrCellsNew [fgkNEMCalCells]; ///<  Array with induced cell energy in T-Card neighbour cells, that before had no signal
  
  Float_t               fTCardCorrInduceEner         [4 ][22]; ///< Induced energy loss gauss constant on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells, param 0  
  Float_t               fTCardCorrInduceEnerFrac     [4 ][22]; ///< Induced energy loss gauss fraction param0 on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells, param 0  
  Float_t               fTCardCorrInduceEnerFracP1   [4 ][22]; ///< Induced energy loss gauss fraction param1 on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells, param1  
  Float_t               fTCardCorrInduceEnerFracWidth[4 ][22]; ///< Induced energy loss gauss witdth on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells  
  Float_t               fTCardCorrInduceEnerFracMax[22];   ///< In case fTCardCorrInduceEnerFracP1  is non null, restrict the maximum fraction of induced energy per SM  
  Float_t               fTCardCorrInduceEnerFracMin[22];   ///< In case fTCardCorrInduceEnerFracP1  is non null, restrict the minimum fraction of induced energy per SM  
  Float_t               fTCardCorrInduceEnerProb[22];      ///< Probability to induce energy loss per SM   
 
 
  TRandom3              fRandom   ;                ///<  Random generator
  Bool_t                fRandomizeTCard ;          ///<  Use random induced energy
  
  Float_t               fTCardCorrMinAmp;          ///<  Minimum cell energy to induce signal on adjacent cells
  Float_t               fTCardCorrMinInduced;      ///<  Minimum induced energy signal on adjacent cells, sum of induced plus original energy, use same as cell energy clusterization cut
  Float_t               fTCardCorrMaxInducedLowE;  ///<  Maximum value of induced energy signal that is always accepted, order of ADC, tipically 10 MeV
  Float_t               fTCardCorrMaxInduced;      ///<  Maximum induced energy signal on adjacent cells
  
  Bool_t                fPrintOnce;                ///< Print once analysis parameters

  Bool_t                fDoMergedBCs;              ///< flag whether to load four histos for the time calib or one merged histo
  Bool_t                fLoad1DRecalibFactors;     ///< Flag to load 1D energy recalibration factors
  
  /// Copy constructor not implemented.
  AliAnalysisTaskEMCALClusterize(           const AliAnalysisTaskEMCALClusterize&) ;
    
  /// Assignment operator not implemented.
  AliAnalysisTaskEMCALClusterize& operator=(const AliAnalysisTaskEMCALClusterize&) ;

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEMCALClusterize, 44) ;
  /// \endcond

};

#endif //ALIANALYSISTASKEMCALCLUSTERIZE_H
