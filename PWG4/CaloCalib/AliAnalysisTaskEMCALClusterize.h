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

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALClusterize : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALClusterize();
  AliAnalysisTaskEMCALClusterize(const char *name);
  virtual ~AliAnalysisTaskEMCALClusterize();
  
 private:
  AliAnalysisTaskEMCALClusterize(const AliAnalysisTaskEMCALClusterize&); // not implemented
  AliAnalysisTaskEMCALClusterize& operator=(const AliAnalysisTaskEMCALClusterize&); // not implemented
  
 public:
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual Bool_t UserNotify();

  void           SetOCDBPath(const char *path)                  { fOCDBpath = path             ; }
  
  //Geometry methods
  void           SetGeometryName(TString &name)                 { fGeomName = name             ; }
  TString        GeometryName() const                           { return fGeomName             ; }  
  void           SwitchOnLoadOwnGeometryMatrices()              { fLoadGeomMatrices = kTRUE    ; }
  void           SwitchOffLoadOwnGeometryMatrices()             { fLoadGeomMatrices = kFALSE   ; }
  void           SetGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fGeomMatrix[i]    = m        ; }

  //AOD methods
  void           SetAODBranchName(TString &name)                { fOutputAODBranchName = name  ; }
  void           FillAODFile(Bool_t yesno)                      { fFillAODFile         = yesno ; }
  
  //Algorithms settings
  void           JustUnfold(Bool_t yesno)                       { fJustUnfold          = yesno ; }
  AliEMCALRecParam * GetRecParam()    const                     { return fRecParam             ; }
  void           InitClusterization();
  
 private:
    
  virtual void  RecPoints2Clusters(TClonesArray *fdigitsArr, TObjArray *fRecPoints, TObjArray *clusArray);
  
  //Geometry  
  AliEMCALGeometry      *fGeom;             //! emcal geometry
  TString                fGeomName;         // Name of geometry to use.
  TGeoHMatrix           *fGeomMatrix[10];   //! Geometry matrices with alignments
  Bool_t                 fGeomMatrixSet;    // Set geometry matrices only once, for the first event.         
  Bool_t                 fLoadGeomMatrices; // Matrices set from configuration, not get from geometry.root or from ESDs/AODs

  //OCDB
  AliEMCALCalibData     *fCalibData;        //! emcal calib data
  AliCaloCalibPedestal  *fPedestalData;     //! emcal pedestal
  TString                fOCDBpath;         // Path with OCDB location

  //Temporal arrays
  TClonesArray          *fDigitsArr;        //-> digits array
  TObjArray             *fClusterArr;       //! recpoints array
  TObjArray             *fCaloClusterArr;   //! CaloClusters array

  //Clusterizers 
  AliEMCALRecParam      *fRecParam;         //! reconstruction parameters container
  AliEMCALClusterizer   *fClusterizer;      //! emcal clusterizer
  AliEMCALAfterBurnerUF *fUnfolder;         //! unfolding procedure
  Bool_t                 fJustUnfold;       // Just unfold, do not recluster
  
  //AOD
  TClonesArray          *fOutputAODBranch;  //-> AOD Branch with output clusters  
  TString                fOutputAODBranchName;  // New of output AOD branch
  Bool_t                 fFillAODFile;      // Fill the output AOD file with the new clusters, 
                                            // if not they will be only available for the event they were generated

  ClassDef(AliAnalysisTaskEMCALClusterize, 1);
};

#endif //ALIANALYSISTASKEMCALCLUSTERIZE_H
