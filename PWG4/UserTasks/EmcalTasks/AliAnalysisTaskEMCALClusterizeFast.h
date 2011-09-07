#ifndef ALIANALYSISTASKEMCALCLUSTERIZEFAST_H
#define ALIANALYSISTASKEMCALCLUSTERIZEFAST_H

// $Id$

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

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALClusterizeFast : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALClusterizeFast();
  AliAnalysisTaskEMCALClusterizeFast(const char *name);
  virtual ~AliAnalysisTaskEMCALClusterizeFast();
  
 public:
  virtual void           UserCreateOutputObjects();
  virtual void           UserExec(Option_t *option);

  Bool_t                 GetAttachClusters()                          const   { return fAttachClusters               ; }
  Bool_t                 GetRecalibrateOnly()                         const   { return fRecalibOnly                  ; }
  Bool_t                 GetSubBackground()                           const   { return fSubBackground                ; }
  const TObjArray       *GetClusters()                                const   { return fClusterArr                   ; }
  const TClonesArray    *GetDigits()                                  const   { return fDigitsArr                    ; }
  const TString         &GeometryName()                               const   { return fGeomName                     ; }  
  AliEMCALRecParam      *GetRecParam()                                const   { return fRecParam                     ; }
  AliEMCALRecoUtils     *GetRecoUtils()                               const   { return fRecoUtils                    ; }
  AliEMCALCalibData     *GetCalibData()                               const   { return fCalibData                    ; }
  AliCaloCalibPedestal  *GetPedData()                                 const   { return fPedestalData                 ; }
  TGeoHMatrix           *GetGeometryMatrix(Int_t i)                   const   { return fGeomMatrix[i]                ; }
  Bool_t                 GetCreatePattern()                           const   { return fCreatePattern                ; }
  Bool_t                 GetOverwrite()                               const   { return fOverwrite                    ; }
  const TString         &GetNewClusterArrayName()                     const   { return fNewClusterArrayName          ; }
  Int_t                  GetnPhi()                                    const   { return fNPhi                         ; }
  Int_t                  GetnEta()                                    const   { return fNEta                         ; }
  Int_t                  GetShiftPhi()                                const   { return fShiftPhi                     ; }
  Int_t                  GetShiftEta()                                const   { return fShiftEta                     ; }
  Bool_t                 GetTRUShift()                                const   { return fTRUShift                     ; }
  Bool_t                 GetClusterizeFastORs()                       const   { return fClusterizeFastORs            ; }
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
  void                   SetRecalibrateCellsOnly(Bool_t b)                    { fRecalibOnly                 = b     ; }
  void                   SetSubBackground(Bool_t b)                           { fSubBackground               = b     ; }
  void                   SetCreatePattern(Bool_t yes)                         { fCreatePattern               = yes   ; if (yes) fOverwrite = kTRUE; }
  void                   SetOverwrite(Bool_t yes)                             { fOverwrite                   = yes   ; if (yes) fOverwrite = kTRUE; }
  void                   SetNewClusterArrayName(TString name)                 { fNewClusterArrayName         = name  ; }
  void                   SetnPhi(Int_t n)                                     { fNPhi                        = n     ; }
  void                   SetnEta(Int_t n)                                     { fNEta                        = n     ; }
  void                   SetShiftPhi(Int_t n)                                 { fShiftPhi                    = n     ; }
  void                   SetShiftEta(Int_t n)                                 { fShiftEta                    = n     ; }
  void                   SetTRUShift(Bool_t yes)                              { fTRUShift                    = yes   ; }
  void                   SetClusterizeFastORs(Bool_t yes)                     { fClusterizeFastORs           = yes   ; if (yes) fOverwrite = kFALSE; }

 protected:
  virtual void           Clusterize();
  virtual void           FillDigitsArray();
  virtual void           Init();
  virtual void           RecPoints2Clusters(TClonesArray *clus);
  virtual void           UpdateCells();
  virtual void           UpdateClusters();

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
  TGeoHMatrix           *fGeomMatrix[12];                 // geometry matrices with alignments
  TString                fOCDBpath;                       // path with OCDB location
  AliEMCALCalibData     *fCalibData;                      // EMCAL calib data
  AliCaloCalibPedestal  *fPedestalData;                   // EMCAL pedestal
  TClonesArray          *fOutputAODBranch;                //!AOD Branch with output clusters  
  TString                fOutputAODBrName;                // output AOD branch name (none by default)
  AliEMCALRecoUtils     *fRecoUtils;                      // access to factorized reconstruction algorithms
  Bool_t                 fLoadCalib;                      // access calib object from OCDB (def=off)
  Bool_t                 fLoadPed;                        // access ped object from OCDB (def=off)
  Bool_t                 fAttachClusters;                 // attach clusters to input event (AOD or ESD)
  Bool_t                 fRecalibOnly;                    // only recalibrate cells if true (def=off)
  Bool_t                 fSubBackground;                  // subtract background if true (def=off)
  Bool_t                 fCreatePattern;                  // removes all cells and creates a cell pattern before running the clusterizer (for debug purposes)
  Bool_t                 fOverwrite;                      // Overwrite existing clusters
  TString                fNewClusterArrayName;            // If not overwriting, name of the new cluster array
  Int_t                  fNPhi;                           // nPhi (for FixedWindowsClusterizer)
  Int_t                  fNEta;                           // nEta (for FixedWinoswsClusterizer)
  Int_t                  fShiftPhi;                       // ShiftPhi (for FixedWindowsClusterizer)
  Int_t                  fShiftEta;                       // ShiftEta (for FixedWindowsClusterizer)
  Bool_t                 fTRUShift;                       // Shifting inside a TRU (true) or through the whole calorimeter (false) (for FixedWindowsClusterizer)
  Bool_t                 fClusterizeFastORs;              // If true, clusterize FastORs instead of cells

 private:
  AliAnalysisTaskEMCALClusterizeFast(const AliAnalysisTaskEMCALClusterizeFast&);            // not implemented
  AliAnalysisTaskEMCALClusterizeFast &operator=(const AliAnalysisTaskEMCALClusterizeFast&); // not implemented

  ClassDef(AliAnalysisTaskEMCALClusterizeFast, 5);
};
#endif //ALIANALYSISTASKEMCALCLUSTERIZEFAST_H
