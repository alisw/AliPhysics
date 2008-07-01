//*********************************************************************
// - AliEVE implementation -
// Fill containers for visualisation of EMCAL data structures
//    - read and store MC Hits
//    - read and store digits from esds or runloader
//    - read and store clusters from esds or runloader 
//
// Author: Magali Estienne (magali.estienne@cern.ch)
// June 30 2008
//*********************************************************************

#ifndef AliEveEMCALData_H
#define AliEveEMCALData_H

#include <TEveUtil.h>
#include <TEveQuadSet.h>
#include <TEveBoxSet.h>
#include <TEvePointSet.h>
#include <TGedFrame.h>
#include <TGeoNode.h>
#include <TGeoMatrix.h>
#include <TClonesArray.h>
#include <TTree.h>

#include <map>
#include <vector>

class AliRunLoader;
class AliESDEvent;
class AliEMCAL;
class AliEMCALLoader;
class AliEMCALGeometry;
class AliEveEMCALSModule;
class AliEveEMCALSModuleData;

class AliEveEMCALData : public TObject, public TEveRefCnt
{
 public:
  AliEveEMCALData();
  AliEveEMCALData(AliRunLoader* rl, TGeoNode* node, TGeoHMatrix* m);
  ~AliEveEMCALData();

  void Reset();
  void SetTree(TTree* tree);
  void SetESD(AliESDEvent* esd);
  void SetNode(TGeoNode* node);
  void InitEMCALGeom(AliRunLoader* rl);
  void GetGeomInfo(Int_t id, Int_t &iSupMod, Double_t& x, Double_t& y, Double_t& z);

  void CreateAllSModules();
  void CreateSModule(Int_t sm);
  void DropAllSModules();
  void DeleteSuperModules();

  void LoadHits(TTree* t);
  void LoadDigits(TTree* t);
  void LoadRecPoints(TTree* t);
  void LoadHitsFromEMCALLoader(AliEMCALLoader* emcl);
  void LoadDigitsFromEMCALLoader(AliEMCALLoader* emcl);
  void LoadRecPointsFromEMCALLoader(AliEMCALLoader* emcl);
  void LoadDigitsFromESD();
  void LoadRecPointsFromESD();
  void LoadRaw();

  AliEveEMCALSModuleData* GetSModuleData(Int_t sm);
  TEvePointSet*           GetPointSetData() {return fPoint;};

 protected:
  AliEMCAL*         fEmcal;     // EMCal data member
  AliEMCALGeometry* fGeom;      // Data member to set/call EMCAL geometry
  TGeoNode*         fNode;      // node
  TGeoHMatrix*      fHMatrix;   // matrix for local to global transformation
  TTree*            fTree;      // tree
  AliESDEvent*      fESD;       // esd
  Int_t             fNsm;       // Total number of Super Modules
  Int_t             fNsmfull;   // Number of full size Super Modules
  Int_t             fNsmhalf;   // Number of half size Super Modules
  std::vector<AliEveEMCALSModuleData*>   fSM;       // vector of fNsm SModules
  std::vector<AliEveEMCALSModuleData*>   fSMfull;   // vector of fNsmfull SModules
  std::vector<AliEveEMCALSModuleData*>   fSMhalf;   // vector of fNhalf SModules
  AliRunLoader*     fRunLoader; // Run Loader
  Int_t             fDebug;     // Debug option
  TEvePointSet*     fPoint;    // TEvePointSet for hits 

 private:
  AliEveEMCALData(const AliEveEMCALData &edata);            
  AliEveEMCALData& operator=(const AliEveEMCALData&); // Not implemented

  ClassDef(AliEveEMCALData, 0); // Base class for TRD hits visualisation
};

#endif
