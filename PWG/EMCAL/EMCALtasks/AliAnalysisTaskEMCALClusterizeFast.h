#ifndef ALIANALYSISTASKEMCALCLUSTERIZEFAST_H
#define ALIANALYSISTASKEMCALCLUSTERIZEFAST_H

//_________________________________________________________________________
///
/// \class AliAnalysisTaskEMCALClusterizeFast
/// \brief EMCal data reclusterization
///
/// This analysis provides a new list of clusters to be used in other analysis
/// Adapted from analysis class AliAnalysisTaksEMCALClusterizeFast from Deepa Thomas
///
/// \author Constantin Loizides, LBNL
/// \author Salvatore Aiola, LBNL
///
//_________________________________________________________________________

class TObjArray;
class TClonesArray;
class AliAODEvent;
class AliESDEvent;
class AliEMCALCalibData;
class AliCaloCalibPedestal;
class AliEMCALClusterizer;
class AliEMCALAfterBurnerUF;
class AliEMCALRecParam;
class AliEMCALRecoUtils;
class AliVCaloCells;
class AliEMCALGeometry;

#include "AliEMCALGeoParams.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALClusterizeFast : public AliAnalysisTaskSE {
 public:
  enum InputCellType {
    kFEEData = 0,
    kFEEDataMCOnly,
    kFEEDataExcludeMC,
    kPattern,
    kL0FastORs,
    kL0FastORsTC,
    kL1FastORs
  };
  

  AliAnalysisTaskEMCALClusterizeFast();
  AliAnalysisTaskEMCALClusterizeFast(const char *name);
  virtual ~AliAnalysisTaskEMCALClusterizeFast();
  
 public:
  virtual void           UserCreateOutputObjects();
  virtual void           UserExec(Option_t *option);

  Bool_t                 GetAttachClusters()                          const   { return fAttachClusters               ; }
  Bool_t                 GetSubBackground()                           const   { return fSubBackground                ; }
  const TObjArray       *GetClusters()                                const   { return fClusterArr                   ; }
  const TClonesArray    *GetDigits()                                  const   { return fDigitsArr                    ; }
  const TString         &GeometryName()                               const   { return fGeomName                     ; }  
  AliEMCALRecParam      *GetRecParam()                                const   { return fRecParam                     ; }
  AliEMCALRecoUtils     *GetRecoUtils()                               const   { return fRecoUtils                    ; }
  AliEMCALCalibData     *GetCalibData()                               const   { return fCalibData                    ; }
  AliCaloCalibPedestal  *GetPedData()                                 const   { return fPedestalData                 ; }
  TGeoHMatrix           *GetGeometryMatrix(Int_t i)                   const   { return fGeomMatrix[i]                ; }
  const TString         &GetCaloClustersName()                        const   { return fCaloClustersName             ; }
  Int_t                  GetnPhi()                                    const   { return fNPhi                         ; }
  Int_t                  GetnEta()                                    const   { return fNEta                         ; }
  Int_t                  GetShiftPhi()                                const   { return fShiftPhi                     ; }
  Int_t                  GetShiftEta()                                const   { return fShiftEta                     ; }
  Bool_t                 GetTRUShift()                                const   { return fTRUShift                     ; }
  InputCellType          GetInputCellType()                           const   { return fInputCellType                ; }
  void                   JustUnfold(Bool_t yesno)                             { fJustUnfold                  = yesno ; }
  void                   LoadOwnGeometryMatrices(Bool_t b)                    { fLoadGeomMatrices            = b     ; }
  void                   SetAODBranchName(const char *name)                   { fOutputAODBrName             = name  ; }
  void                   SetAttachClusters(Bool_t b)                          { fAttachClusters              = b     ; }
  void                   SetCalibData(AliEMCALCalibData *d)                   { fCalibData                   = d     ; }
  void                   SetEMCALRecoUtils(AliEMCALRecoUtils *ru)             { fRecoUtils                   = ru    ; }
  void                   SetGeometryMatrix(TGeoHMatrix* m, Int_t i)           { fGeomMatrix[i]               = m     ; }
  void                   SetGeometryName(const char *name)                    { fGeomName                    = name  ; }
  void                   SetLoadCalib(Bool_t b)                               { fLoadCalib                   = b     ; }
  void                   SetLoadPed(Bool_t b)                                 { fLoadPed                     = b     ; }
  void                   SetOCDBPath(const char *path)                        { fOCDBpath                    = path  ; }
  void                   SetPedestalData(AliCaloCalibPedestal *d)             { fPedestalData                = d     ; }
  void                   SetSubBackground(Bool_t b)                           { fSubBackground               = b     ; }
  void                   SetnPhi(Int_t n)                                     { fNPhi                        = n     ; }
  void                   SetnEta(Int_t n)                                     { fNEta                        = n     ; }
  void                   SetShiftPhi(Int_t n)                                 { fShiftPhi                    = n     ; }
  void                   SetShiftEta(Int_t n)                                 { fShiftEta                    = n     ; }
  void                   SetTRUShift(Bool_t yes)                              { fTRUShift                    = yes   ; }
  void                   SetInputCellType(InputCellType ic)                   { fInputCellType               = ic    ; }
  void                   SetTrackName(const char *n)                          { fTrackName                   = n     ; }
  void                   SetCaloClustersName(const char *name)                { fCaloClustersName            = name  ; }
  void                   SetCaloCellsName(const char *name)                   { fCaloCellsName               = name  ; }
  void                   SetUpdateCells(Bool_t b)                             { fDoUpdateCells               = b     ; }
  void                   SetClusterize(Bool_t b)                              { fDoClusterize                = b     ; }
  void                   SetClusterBadChannelCheck(Bool_t b)                  { fClusterBadChannelCheck      = b     ; }
  void                   SetRejectExoticClusters(Bool_t b)                    { fRejectExoticClusters        = b     ; }
  void                   SetRejectExoticCells(Bool_t b)                       { fRejectExoticCells           = b     ; }
  void                   SetFiducial(Bool_t b)                                { fFiducial                    = b     ; }
  void                   SetDoNonLinearity(Bool_t b)                          { fDoNonLinearity              = b     ; }
  void                   SetRecalDistToBadChannels(Bool_t b)                  { fRecalDistToBadChannels      = b     ; }
  void                   SetCellMCLabelFromCluster(Int_t s)                   { fSetCellMCLabelFromCluster   = s     ; }

  // For backward compatibility
  const TString         &GetNewClusterArrayName()                     const   { return GetCaloClustersName()         ; }
  void                   SetNewClusterArrayName(const char *name)             { SetCaloClustersName(name)            ; }
  void                   SetOverwrite(Bool_t b)                               { if (b) SetCaloClustersName(""); else SetCaloClustersName("newCaloClusters");}
  void                   SetRecalibrateCellsOnly(Bool_t b)                    { if (b) { SetUpdateCells(kTRUE); SetClusterize(kFALSE);} else { SetClusterize(kTRUE); } }
  Bool_t                 GetRecalibrateOnly()                         const   { return (Bool_t)(fDoUpdateCells && !fDoClusterize); }
  Bool_t                 GetOverwrite()                               const   { return fCaloClustersName.IsNull()                ; }
    
 protected:
  Bool_t                 AcceptCell(Int_t cellNumber);
  virtual void           Clusterize();
  virtual void           FillDigitsArray();
  virtual void           Init();
  virtual void           RecPoints2Clusters(TClonesArray *clus);
  virtual void           UpdateCells();
  virtual void           UpdateClusters();
  virtual void           CalibrateClusters();
  virtual void           TrackClusterMatching(AliVCluster *c, TClonesArray *tarr);
  virtual void           CopyClusters(TClonesArray *orig, TClonesArray *dest);

  Int_t                  fRun;                            //!run number
  TClonesArray          *fDigitsArr;                      //!digits array
  TObjArray             *fClusterArr;                     //!recpoints array
  AliEMCALRecParam      *fRecParam;                       // reconstruction parameters container
  AliEMCALClusterizer   *fClusterizer;                    //!clusterizer
  AliEMCALAfterBurnerUF *fUnfolder;                       //!unfolding procedure
  Bool_t                 fJustUnfold;                     // just unfold, do not recluster
  TString                fGeomName;                       // name of geometry to use.
  Bool_t                 fGeomMatrixSet;                  // set geometry matrices only once, for the first event.         
  Bool_t                 fLoadGeomMatrices;               // matrices from configuration, not geometry.root nor ESDs/AODs
  TGeoHMatrix           *fGeomMatrix[AliEMCALGeoParams::fgkEMCALModules]; // geometry matrices with alignments
  TString                fOCDBpath;                       // path with OCDB location
  AliEMCALCalibData     *fCalibData;                      // EMCAL calib data
  AliCaloCalibPedestal  *fPedestalData;                   // EMCAL pedestal
  TClonesArray          *fOutputAODBranch;                //!AOD Branch with output clusters  
  TString                fOutputAODBrName;                // output AOD branch name (none by default)
  AliEMCALRecoUtils     *fRecoUtils;                      // access to factorized reconstruction algorithms
  Bool_t                 fLoadCalib;                      // access calib object from OCDB (def=off)
  Bool_t                 fLoadPed;                        // access ped object from OCDB (def=off)
  Bool_t                 fAttachClusters;                 // attach clusters to input event (AOD or ESD)
  Bool_t                 fSubBackground;                  // subtract background if true (def=off)
  Int_t                  fNPhi;                           // nPhi (for FixedWindowsClusterizer)
  Int_t                  fNEta;                           // nEta (for FixedWinoswsClusterizer)
  Int_t                  fShiftPhi;                       // shift in phi (for FixedWindowsClusterizer)
  Int_t                  fShiftEta;                       // shift in eta (for FixedWindowsClusterizer)
  Bool_t                 fTRUShift;                       // shifting inside a TRU (true) or through the whole calorimeter (false) (for FixedWindowsClusterizer)
  InputCellType          fInputCellType;                  // input cells type to make clusters
  TString                fTrackName;                      // if not null use track collection for track/cluster matching
  TString                fCaloCellsName;                  // name of calo cells object
  TString                fCaloClustersName;               // name of calo cluster collection
  Bool_t                 fDoUpdateCells;                  // recalibrate cells
  Bool_t                 fDoClusterize;                   // clusterize
  Bool_t                 fClusterBadChannelCheck;         // cluster bad channel check
  Bool_t                 fRejectExoticClusters;           // reject exotic cluster
  Bool_t                 fRejectExoticCells;              // reject exotic cells on-the-fly
  Bool_t                 fFiducial;                       // fiducial cut
  Bool_t                 fDoNonLinearity;                 // non linearity calib
  Bool_t                 fRecalDistToBadChannels;         // recalculate distance to bad channel
    
  // MC labels
  static const Int_t     fgkTotalCellNumber = 17664 ;     // Maximum number of cells in EMCAL/DCAL: (48*24)*(10+4/3.+6*2/3.)
  
  Int_t                  fOrgClusterCellId[fgkTotalCellNumber]; // Array ID of cluster to wich the cell belongs in unmodified clusters. 
  Int_t                  fCellLabels      [fgkTotalCellNumber]; // Array with MC label/map 

  Int_t                  fSetCellMCLabelFromCluster;      // Use cluster MC label as cell label:
                                                          // 0 - get the MC label stored in cells (not available for productions done with aliroot < v5-02-Rev09)
                                                          // 1 - assign to the cell the MC label of the cluster
  
  Bool_t                 fSetCellMCLabelFromEdepFrac;     // For MC generated with aliroot > v5-07-21, check the EDep information 
                                                          // stored in ESDs/AODs to set the cell MC labels
  
  AliVCaloCells         *fCaloCells;                      //!calo cells object
  TClonesArray          *fCaloClusters;                   //!calo clusters array       
  AliESDEvent           *fEsd;                            //!esd event
  AliAODEvent           *fAod;                            //!aod event
  AliEMCALGeometry      *fGeom;                           //!geometry object

 private:
  AliAnalysisTaskEMCALClusterizeFast(const AliAnalysisTaskEMCALClusterizeFast&);            // not implemented
  AliAnalysisTaskEMCALClusterizeFast &operator=(const AliAnalysisTaskEMCALClusterizeFast&); // not implemented

  ClassDef(AliAnalysisTaskEMCALClusterizeFast, 10);
};
#endif //ALIANALYSISTASKEMCALCLUSTERIZEFAST_H
