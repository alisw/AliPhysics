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

  Bool_t                 GetAttachClusters()                  const  { return fAttachClusters       ; }
  Bool_t                 GetRecalibrateOnly()                 const  { return fRecalibOnly          ; }
  Bool_t                 GetSubBackground()                   const  { return fSubBackground        ; }
  const TObjArray       *GetClusters()                        const  { return fClusterArr           ; }
  const TString         &GeometryName()                       const  { return fGeomName             ; }  
  AliEMCALRecParam      *GetRecParam()                        const  { return fRecParam             ; }
  AliEMCALRecoUtils     *GetRecoUtils()                       const  { return fRecoUtils            ; }
  AliEMCALCalibData     *GetCalibData()                       const  { return fCalibData            ; }
  AliCaloCalibPedestal  *GetPedData()                         const  { return fPedestalData         ; }
  TGeoHMatrix           *GetGeometryMatrix(Int_t i)           const  { return fGeomMatrix[i]        ; }
  void                   JustUnfold(Bool_t yesno)                    { fJustUnfold          = yesno ; }
  void                   LoadOwnGeometryMatrices(Bool_t b)           { fLoadGeomMatrices    = b     ; }
  void                   SetAODBranchName(const char *name)          { fOutputAODBrName     = name  ; }
  void                   SetAttachClusters(Bool_t b)                 { fAttachClusters      = b     ; }
  void                   SetCalibData(AliEMCALCalibData *d)          { fCalibData           = d     ; }
  void                   SetEMCALRecoUtils(AliEMCALRecoUtils *ru)    { fRecoUtils           = ru    ; }
  void                   SetGeometryMatrix(TGeoHMatrix* m, Int_t i)  { fGeomMatrix[i]       = m     ; }
  void                   SetGeometryName(const char *name)           { fGeomName            = name  ; }
  void                   SetLoadCalib(Bool_t b)                      { fLoadCalib           = b     ; }
  void                   SetLoadPed(Bool_t b)                        { fLoadPed             = b     ; }
  void                   SetOCDBPath(const char *path)               { fOCDBpath            = path  ; }
  void                   SetPedestalData(AliCaloCalibPedestal *d)    { fPedestalData        = d     ; }
  void                   SetRecalibrateCellsOnly(Bool_t b)           { fRecalibOnly         = b     ; }
  void                   SetSubBackground(Bool_t b)                  { fSubBackground       = b     ; }

 protected:
  void                   Clusterize();
  void                   FillDigitsArray();
  void                   Init();
  void                   RecPoints2Clusters(TClonesArray *clus);
  void                   SubBackground();
  void                   UpdateCells();
  void                   UpdateClusters();

  Int_t                  fRun;              //!run number
  TClonesArray          *fDigitsArr;        //!digits array
  TObjArray             *fClusterArr;       //!recpoints array
  AliEMCALRecParam      *fRecParam;         // reconstruction parameters container
  AliEMCALClusterizer   *fClusterizer;      //!clusterizer
  AliEMCALAfterBurnerUF *fUnfolder;         //!unfolding procedure
  Bool_t                 fJustUnfold;       // just unfold, do not recluster
  TString                fGeomName;         // name of geometry to use.
  Bool_t                 fGeomMatrixSet;    // set geometry matrices only once, for the first event.         
  Bool_t                 fLoadGeomMatrices; // matrices from configuration, not geometry.root nor ESDs/AODs
  TGeoHMatrix           *fGeomMatrix[12];   // geometry matrices with alignments
  TString                fOCDBpath;         // path with OCDB location
  AliEMCALCalibData     *fCalibData;        // EMCAL calib data
  AliCaloCalibPedestal  *fPedestalData;     // EMCAL pedestal
  TClonesArray          *fOutputAODBranch;  //!AOD Branch with output clusters  
  TString                fOutputAODBrName;  // output AOD branch name (none by default)
  AliEMCALRecoUtils     *fRecoUtils;        // access to factorized reconstruction algorithms
  Bool_t                 fLoadCalib;        // access calib object from OCDB (def=off)
  Bool_t                 fLoadPed;          // access ped object from OCDB (def=off)
  Bool_t                 fAttachClusters;   // attach clusters to input event (AOD or ESD)
  Bool_t                 fRecalibOnly;      // only recalibrate cells if true (def=off)
  Bool_t                 fSubBackground;    // subtract background if true (def=off)

 private:
  AliAnalysisTaskEMCALClusterizeFast(const AliAnalysisTaskEMCALClusterizeFast&);            // not implemented
  AliAnalysisTaskEMCALClusterizeFast &operator=(const AliAnalysisTaskEMCALClusterizeFast&); // not implemented

  ClassDef(AliAnalysisTaskEMCALClusterizeFast, 5);
};
#endif //ALIANALYSISTASKEMCALCLUSTERIZEFAST_H
