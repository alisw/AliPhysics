#ifndef ALIANALYSISTASKEMCALCLUSTERIZE_H
#define ALIANALYSISTASKEMCALCLUSTERIZE_H

// This analysis provides a new list of clusters to be used in other analysis
// Author: Gustavo Conesa Balbastre,
//         Adapted from analysis class from Deepa Thomas

//Root
class TTree;
class TClonesArray;

#include "AliCentrality.h"

//EMCAL
class AliEMCALGeometry;
class AliEMCALCalibData;
class AliCaloCalibPedestal;
class AliEMCALClusterizer;
class AliEMCALAfterBurnerUF;
class AliEMCALRecPoint;
class AliAODCaloCluster;
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
  virtual void   LocalInit()                                    { Init()                       ; }
    
  // Event methods, settings
  
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

  void           SetImportGeometryFromFile(Bool_t  im, 
                                           TString pa = "")     { fImportGeometryFromFile = im ; 
                                                                  fImportGeometryFilePath = pa ; }    
  //AOD methods
  void           SetAODBranchName(TString &name)                { fOutputAODBranchName = name  ; }
  void           FillAODFile(Bool_t yesno)                      { fFillAODFile         = yesno ; }
  void           FillAODCaloCells();
  void           FillAODHeader();
  void           SwitchOnFillAODHeader()                        { fFillAODHeader     = kTRUE   ; }
  void           SwitchOffFillAODHeader()                       { fFillAODHeader     = kFALSE  ; } 
  void           SwitchOnFillAODCaloCells()                     { fFillAODCaloCells  = kTRUE   ; }
  void           SwitchOffFillAODCaloCells()                    { fFillAODCaloCells  = kFALSE  ; } 

  void           SwitchOnRecalibrateWithClusterTime()           { fRecalibrateWithClusterTime  = kTRUE   ; }
  void           SwitchOffRecalibrateWithClusterTime()          { fRecalibrateWithClusterTime  = kFALSE  ; }
  

  
  //Algorithms settings
  
  AliEMCALRecParam * GetRecParam()                              { if(!fRecParam)  fRecParam  = new AliEMCALRecParam  ;
                                                                  return fRecParam             ; }
  
  AliEMCALRecoUtils* GetRecoUtils()                             { if(!fRecoUtils) fRecoUtils = new AliEMCALRecoUtils ;  
                                                                  return fRecoUtils            ; }

  void           InitClusterization();
  void           ClusterizeCells();
  void           ClusterUnfolding();
  void           JustUnfold(Bool_t yesno)                       { fJustUnfold        = yesno   ; }
    
  void           SetConfigFileName(TString name)                { fConfigName        = name    ; }
  void           SetMaxEvent(Int_t max)                         { fMaxEvent          = max     ; }
  
  void           SwitchOnTrackMatching()                        { fDoTrackMatching   = kTRUE   ; }
  void           SwitchOffTrackMatching()                       { fDoTrackMatching   = kFALSE  ; } 

  // Cell selection after unfolding
  void           SwitchOnCellEnergySelection()                  { fSelectCell        = kTRUE   ; }
  void           SwitchOffCellEnergySelection()                 { fSelectCell        = kFALSE  ; } 
  void           SetCellCuts(Float_t e, Float_t frac)           { fSelectCellMinE    = e       ; 
                                                                  fSelectCellMinFrac = frac    ; }  
  // OADB options settings
  
  void           AccessOADB() ;
  
  TString        GetPass()    ;
  
  void           SwitchOnEMCALOADB()                            { fAccessOADB        = kTRUE   ; }
  void           SwitchOffEMCALOADB()                           { fAccessOADB        = kFALSE  ; }
    
  void           SetOADBFilePath(TString path)                  { fOADBFilePath      = path    ; }
  
  // Centrality selection
  
  AliCentrality* GetCentrality()                                { return InputEvent()->GetCentrality() ; } //Look in AOD reader, different there
  void           SetCentralityClass(TString name)               { fCentralityClass   = name            ; }
  TString        GetCentralityClass()                     const { return fCentralityClass              ; }
  Float_t        GetEventCentrality()                           { if(GetCentrality()) return GetCentrality()->GetCentralityPercentile(fCentralityClass) ;
                                                                  else                return -1.       ; }
  void           SetCentralityBin(Int_t min, Int_t max) //Set the centrality bin to select the event. If used, then need to get percentile
                                                                { fCentralityBin[0]=min ; fCentralityBin[1]=max ; }
  Float_t        GetCentralityBin(Int_t i)                const { if(i < 0 || i > 1) return -1 ; 
                                                                  else               return fCentralityBin[i]   ; }
  
  //MC label methods
  
  void           RemapMCLabelForAODs(Int_t &label);
  void           SwitchOnRemapMCLabelForAODs()                  { fRemapMCLabelForAODs  = kTRUE   ; }
  void           SwitchOffRemapMCLabelForAODs()                 { fRemapMCLabelForAODs  = kFALSE  ; }

  void           SetClustersMCLabelFrom2SelectedLabels(AliEMCALRecPoint* recPoint, AliAODCaloCluster *clus) ;
  void           SetClustersMCLabelFromOriginalClusters(AliAODCaloCluster * clus) ;
  
  void           SwitchOnUseClusterMCLabelForCell(Int_t opt = 2) { fSetCellMCLabelFromCluster = opt ; }
  void           SwitchOffUseClusterMCLabelForCell()             { fSetCellMCLabelFromCluster = 0   ; }

private:
    
  virtual void   FillCaloClusterInEvent();
  
  virtual void   RecPoints2Clusters();
  
  virtual void   ResetArrays();
  
  AliVEvent             *fEvent;                   // Event 
  
  //Geometry  
  AliEMCALGeometry      *fGeom;                    // EMCAL geometry
  TString                fGeomName;                // Name of geometry to use.
  TGeoHMatrix           *fGeomMatrix[12];          // Geometry matrices with alignments
  Bool_t                 fGeomMatrixSet;           // Set geometry matrices only once, for the first event.         
  Bool_t                 fLoadGeomMatrices;        // Matrices set from configuration, not get from geometry.root or from ESDs/AODs

  //OCDB
  AliEMCALCalibData     *fCalibData;               // EMCAL calib data
  AliCaloCalibPedestal  *fPedestalData;            // EMCAL pedestal
  TString                fOCDBpath;                // Path with OCDB location
  Bool_t                 fAccessOCDB;              // Need to access info from OCDB (not really)   

  //Temporal arrays
  TClonesArray          *fDigitsArr;               //! Digits array
  TObjArray             *fClusterArr;              //! Recpoints array
  TObjArray             *fCaloClusterArr;          //! CaloClusters array

  //Clusterizers 
  AliEMCALRecParam      *fRecParam;                // Reconstruction parameters container
  AliEMCALClusterizer   *fClusterizer;             //! EMCAL clusterizer
  AliEMCALAfterBurnerUF *fUnfolder;                //! Unfolding procedure
  Bool_t                 fJustUnfold;              // Just unfold, do not recluster
  
  //AOD
  TClonesArray          *fOutputAODBranch;         //! AOD Branch with output clusters  
  TString                fOutputAODBranchName;     // New of output AOD branch
  Bool_t                 fOutputAODBranchSet ;     // Set the AOD clusters branch in the input event once
  Bool_t                 fFillAODFile;             // Fill the output AOD file with the new clusters, 
                                                   // if not they will be only available for the event they were generated
  Bool_t                 fFillAODHeader;           // Copy header to standard branch
  Bool_t                 fFillAODCaloCells;        // Copy calocells to standard branch

  Int_t                  fRun;                     // run number
  
  AliEMCALRecoUtils*     fRecoUtils;               // Access to factorized reconstruction algorithms
  TString                fConfigName;              // Name of analysis configuration file
  
  
  Int_t                  fOrgClusterCellId[12672]; // Array ID of cluster to wich the cell belongs in unmodified clusters
  Int_t                  fCellLabels[12672];       // Array with MC label to be passed to digit.
  Int_t                  fCellSecondLabels[12672]; // Array with Second MC label to be passed to digit. 
  Double_t               fCellTime[12672];         // Array with cluster time to be passed to digit in case of AODs 
  Float_t                fCellMatchdEta[12672];    // Array with cluster-track dPhi 
  Float_t                fCellMatchdPhi[12672];    // Array with cluster-track dEta 

  Bool_t                 fRecalibrateWithClusterTime; // Use fCellTime to store time of cells in cluster
  
  Int_t                  fMaxEvent;                // Set a maximum event
  
  Bool_t                 fDoTrackMatching;         // On/Off the matching recalulation to speed up analysis in PbPb
  Bool_t                 fSelectCell;              // Reject cells from cluster if energy is too low and recalculate position/energy and other
  Float_t                fSelectCellMinE;          // Min energy cell threshold, after unfolding
  Float_t                fSelectCellMinFrac;       // Min fraction of cell energy after unfolding cut
  Bool_t                 fRemoveLEDEvents;         // Remove LED events, use only for LHC11a 
  Bool_t                 fRemoveExoticEvents;      // Remove exotic events
  
  Bool_t                 fImportGeometryFromFile;  // Import geometry settings in geometry.root file
  TString                fImportGeometryFilePath;  // path fo geometry.root file

  Bool_t                 fOADBSet ;                // AODB parameters already set
  Bool_t                 fAccessOADB ;             // Get calibration from OADB for EMCAL
  TString                fOADBFilePath ;           // Default path $ALICE_ROOT/OADB/EMCAL, if needed change
    
  //Centrality
  TString                fCentralityClass;         // Name of selected centrality class     
  Float_t                fCentralityBin[2];        // Minimum and maximum value of the centrality for the analysis
  
  // Event selection with some signal in EMCAL
  Bool_t                 fSelectEMCALEvent;       //  Process the event if there is some high energy cluster 
  Float_t                fEMCALEnergyCut;         //  At least an EMCAL cluster with this energy in the event
  Int_t                  fEMCALNcellsCut;         //  At least an EMCAL cluster with fNCellsCut cells over fEnergyCut

  Int_t                  fSetCellMCLabelFromCluster; // Use cluster MC label as cell label:
                                                     // 0 - get the MC label stored in cells
                                                     // 1 - from old way, select 2 most likely labels
                                                     // 2 - from new way, get the original clusters, add all the MC labels (useful for any reclusterization with output V1 clusters)
  Bool_t                 fRemapMCLabelForAODs ;      // Remap AOD cells MC label

  
  Bool_t                 fInputFromFilter ;          // Get the input from AODs from the filter 
  
  AliAnalysisTaskEMCALClusterize(           const AliAnalysisTaskEMCALClusterize&); // not implemented
  AliAnalysisTaskEMCALClusterize& operator=(const AliAnalysisTaskEMCALClusterize&); // not implemented

  ClassDef(AliAnalysisTaskEMCALClusterize, 27);

};

#endif //ALIANALYSISTASKEMCALCLUSTERIZE_H
