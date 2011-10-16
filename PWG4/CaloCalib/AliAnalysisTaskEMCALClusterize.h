#ifndef ALIANALYSISTASKEMCALCLUSTERIZE_H
#define ALIANALYSISTASKEMCALCLUSTERIZE_H

// This analysis provides a new list of clusters to be used in other analysis
// Author: Gustavo Conesa Balbastre,
//         Adapted from analysis class from Deepa Thomas

//Root
class TTree;
class TClonesArray;

//EMCAL
class AliEMCALGeometry;
class AliEMCALCalibData;
class AliCaloCalibPedestal;
class AliEMCALClusterizer;
class AliEMCALAfterBurnerUF;
class AliEMCALRecParam;
class AliEMCALRecoUtils;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALClusterize : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALClusterize();
  AliAnalysisTaskEMCALClusterize(const char *name);
  virtual ~AliAnalysisTaskEMCALClusterize();

 private:  
  AliAnalysisTaskEMCALClusterize(const AliAnalysisTaskEMCALClusterize&); 
  AliAnalysisTaskEMCALClusterize& operator=(const AliAnalysisTaskEMCALClusterize&); // not implemented
  
 public:
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Init();
  virtual void   LocalInit()                                    { Init()                       ; }
    
  // Event methods, settings
  void           CheckAndGetEvent();
  
  Bool_t         IsExoticEvent();
  void           SwitchOnExoticEventsRemoval()                  { fRemoveExoticEvents= kTRUE   ; }
  void           SwitchOffExoticEventsRemoval()                 { fRemoveExoticEvents= kFALSE  ; } 
  
  Bool_t         IsExoticCell(const Int_t absId, const Float_t ecell, 
                              const Float_t tcell, AliVCaloCells* cells);
  void           SwitchOnExoticCellRemoval()                    { fRemoveExoticCells = kTRUE   ; }
  void           SwitchOffExoticCellRemoval()                   { fRemoveExoticCells = kFALSE  ; } 
  
  void           SetExoticCellFractionCut(Float_t f)            { fExoticCellFraction = f      ; }
  void           SetExoticCellDiffTimeCut(Float_t dt)           { fExoticCellDiffTime = dt     ; }
  void           SetExoticCellMinAmplitudeCut(Float_t ma)       { fExoticCellMinAmplitude = ma ; }
  
  Bool_t         IsLEDEvent();
  void           SwitchOnLEDEventsRemoval()                     { fRemoveLEDEvents   = kTRUE   ; }
  void           SwitchOffLEDEventsRemoval()                    { fRemoveLEDEvents   = kFALSE  ; } 

  Bool_t         AcceptCalibrateCell(const Int_t absId, Float_t & amp, Double_t & time, AliVCaloCells* cells) ;
  
  //OCDB
  Bool_t         AccessOCDB();
  void           SwitchOnAccessOCDB()                           { fAccessOCDB       = kTRUE    ; }
  void           SwitchOffAccessOCDB()                          { fAccessOCDB       = kFALSE   ; } 
  void           SetOCDBPath(const char *path)                  { fOCDBpath         = path     ; }
  
  //Geometry methods
  void           InitGeometry();
  void           SetGeometryName(TString &name)                 { fGeomName = name             ; }
  TString        GeometryName()                          const  { return fGeomName             ; }  
  void           SwitchOnLoadOwnGeometryMatrices()              { fLoadGeomMatrices = kTRUE    ; }
  void           SwitchOffLoadOwnGeometryMatrices()             { fLoadGeomMatrices = kFALSE   ; } 
  void           SetGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fGeomMatrix[i]    = m        ; }

  //AOD methods
  void           SetAODBranchName(TString &name)                { fOutputAODBranchName = name  ; }
  void           FillAODFile(Bool_t yesno)                      { fFillAODFile         = yesno ; }
  void           FillAODCaloCells();
  void           FillAODHeader();
  void           SwitchOnFillAODHeader()                        { fFillAODHeader     = kTRUE   ; }
  void           SwitchOffFillAODHeader()                       { fFillAODHeader     = kFALSE  ; } 
  void           SwitchOnFillAODCaloCells()                     { fFillAODCaloCells  = kTRUE   ; }
  void           SwitchOffFillAODCaloCells()                    { fFillAODCaloCells  = kFALSE  ; } 
  
  //Algorithms settings
  AliEMCALRecParam * GetRecParam()                       const  { return fRecParam             ; }
  void           InitClusterization();
  void           ClusterizeCells();
  void           ClusterUnfolding();
  void           JustUnfold(Bool_t yesno)                       { fJustUnfold        = yesno   ; }
  
  void           SetEMCALRecoUtils(AliEMCALRecoUtils * ru)      { fRecoUtils         = ru      ; }
  AliEMCALRecoUtils* GetRecoUtils()                      const  { return fRecoUtils            ; }
  
  void           SetConfigFileName(TString name)                { fConfigName        = name    ; }
  void           SetMaxEvent(Int_t max)                         { fMaxEvent          = max     ; }
  
  void           SwitchOnTrackMatching()                        { fDoTrackMatching   = kTRUE   ; }
  void           SwitchOffTrackMatching()                       { fDoTrackMatching   = kFALSE  ; } 

  // Cell selection after unfolding
  void           SwitchOnCellEnergySelection()                  { fSelectCell        = kTRUE   ; }
  void           SwitchOffCellEnergySelection()                 { fSelectCell        = kFALSE  ; } 
  void           SetCellCuts(Float_t e, Float_t frac)           { fSelectCellMinE    = e       ; 
                                                                  fSelectCellMinFrac = frac    ; }  
  
 private:
    
  virtual void  RecPoints2Clusters(TClonesArray *fdigitsArr, TObjArray *fRecPoints, TObjArray *clusArray);
  
  AliVEvent             *fEvent;                   // Event 
  
  //Geometry  
  AliEMCALGeometry      *fGeom;                    // EMCAL geometry
  TString                fGeomName;                // Name of geometry to use.
  TGeoHMatrix           *fGeomMatrix[10];          // Geometry matrices with alignments
  Bool_t                 fGeomMatrixSet;           // Set geometry matrices only once, for the first event.         
  Bool_t                 fLoadGeomMatrices;        // Matrices set from configuration, not get from geometry.root or from ESDs/AODs

  //OCDB
  AliEMCALCalibData     *fCalibData;               // EMCAL calib data
  AliCaloCalibPedestal  *fPedestalData;            // EMCAL pedestal
  TString                fOCDBpath;                // Path with OCDB location
  Bool_t                 fAccessOCDB;              // Need to access info from OCDB (not really)   

  //Temporal arrays
  TClonesArray          *fDigitsArr;               //-> Digits array
  TObjArray             *fClusterArr;              //-> Recpoints array
  TObjArray             *fCaloClusterArr;          //-> CaloClusters array

  //Clusterizers 
  AliEMCALRecParam      *fRecParam;                // Reconstruction parameters container
  AliEMCALClusterizer   *fClusterizer;             //! EMCAL clusterizer
  AliEMCALAfterBurnerUF *fUnfolder;                //! Unfolding procedure
  Bool_t                 fJustUnfold;              // Just unfold, do not recluster
  
  //AOD
  TClonesArray          *fOutputAODBranch;         //! AOD Branch with output clusters  
  TString                fOutputAODBranchName;     // New of output AOD branch
  Bool_t                 fFillAODFile;             // Fill the output AOD file with the new clusters, 
                                                   // if not they will be only available for the event they were generated
  Bool_t                 fFillAODHeader;           // Copy header to standard branch
  Bool_t                 fFillAODCaloCells;        // Copy calocells to standard branch

  Int_t                  fRun;                     //!run number
  
  AliEMCALRecoUtils*     fRecoUtils;               // Access to factorized reconstruction algorithms
  TString                fConfigName;              // Name of analysis configuration file
  
  Int_t                  fCellLabels[12672];       // Array with MC label to be passed to digit. 
  Int_t                  fCellSecondLabels[12672]; // Array with Second MC label to be passed to digit. 
  Double_t               fCellTime[12672];         // Array with cluster time to be passed to digit in case of AODs 

  Int_t                  fMaxEvent;                // Set a maximum event
  
  Bool_t                 fDoTrackMatching;         // On/Off the matching recalulation to speed up analysis in PbPb
  Bool_t                 fSelectCell;              // Reject cells from cluster if energy is too low and recalculate position/energy and other
  Float_t                fSelectCellMinE;          // Min energy cell threshold, after unfolding
  Float_t                fSelectCellMinFrac;       // Min fraction of cell energy after unfolding cut
  Bool_t                 fRemoveLEDEvents;         // Remove LED events, use only for LHC11a 
  Bool_t                 fRemoveExoticEvents;      // Remove exotic events
  Bool_t                 fRemoveExoticCells;       // Remove exotic cells
  Float_t                fExoticCellFraction;      // Good cell if fraction < 1-ecross/ecell
  Float_t                fExoticCellDiffTime;      // If time of candidate to exotic and close cell is too different, it must be noisy, set amp to 0
  Float_t                fExoticCellMinAmplitude;  // Check for exotic only if amplitud is larger than this value
  
  ClassDef(AliAnalysisTaskEMCALClusterize, 14);

};

#endif //ALIANALYSISTASKEMCALCLUSTERIZE_H
