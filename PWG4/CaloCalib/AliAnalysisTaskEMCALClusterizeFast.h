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
  AliAnalysisTaskEMCALClusterizeFast(const char *name=0);
  virtual ~AliAnalysisTaskEMCALClusterizeFast();
  
 public:
  virtual void           UserCreateOutputObjects();
  virtual void           UserExec(Option_t *option);

  const TObjArray       *GetClusters()                        const  { return fClusterArr           ; }
  const TString         &GeometryName()                       const  { return fGeomName             ; }  
  AliEMCALRecParam      *GetRecParam()                        const  { return fRecParam             ; }
  AliEMCALRecoUtils     *GetRecoUtils()                       const  { return fRecoUtils            ; }
  void                   JustUnfold(Bool_t yesno)                    { fJustUnfold          = yesno ; }
  void                   LoadOwnGeometryMatrices(Bool_t b)           { fLoadGeomMatrices    = b     ; }
  void                   SetAODBranchName(const char *name)          { fOutputAODBrName     = name  ; }
  void                   SetCalibData(AliEMCALCalibData *d)          { fCalibData           = d     ; }
  void                   SetEMCALRecoUtils(AliEMCALRecoUtils *ru)    { fRecoUtils           = ru    ; }
  void                   SetGeometryMatrix(TGeoHMatrix* m, Int_t i)  { fGeomMatrix[i]       = m     ; }
  void                   SetGeometryName(const char *name)           { fGeomName            = name  ; }
  void                   SetOCDBPath(const char *path)               { fOCDBpath            = path  ; }
  void                   SetPedestalData(AliCaloCalibPedestal *d)    { fPedestalData        = d     ; }

 private:
  AliAnalysisTaskEMCALClusterizeFast(const AliAnalysisTaskEMCALClusterizeFast&);            // not implemented
  AliAnalysisTaskEMCALClusterizeFast &operator=(const AliAnalysisTaskEMCALClusterizeFast&); // not implemented

  void                   Init();
  void                   FillDigitsArray(AliAODEvent *event);
  void                   FillDigitsArray(AliESDEvent *event);
  void                   RecPoints2Clusters();

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
  TGeoHMatrix           *fGeomMatrix[10];   // geometry matrices with alignments
  TString                fOCDBpath;         // path with OCDB location
  AliEMCALCalibData     *fCalibData;        // EMCAL calib data
  AliCaloCalibPedestal  *fPedestalData;     // EMCAL pedestal
  TClonesArray          *fOutputAODBranch;  //!AOD Branch with output clusters  
  TString                fOutputAODBrName;  // output AOD branch name (none by default)
  AliEMCALRecoUtils     *fRecoUtils;        // access to factorized reconstruction algorithms
  
  ClassDef(AliAnalysisTaskEMCALClusterizeFast, 1);
};
#endif //ALIANALYSISTASKEMCALCLUSTERIZEFAST_H
