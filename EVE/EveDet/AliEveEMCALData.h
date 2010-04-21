//
// Fill containers for visualisation of EMCAL data structures
//
// Author: Magali Estienne (magali.estienne@cern.ch)
// June 30 2008
//

#ifndef ALIEVEEMCALDATA_H
#define ALIEVEEMCALDATA_H

#include <TGeoNode.h> 
#include <TGeoMatrix.h>
#include <TEveRGBAPalette.h>
#include <TEveTrans.h>
#include <TEveQuadSet.h> 
#include "AliESDEvent.h"
#include "AliRun.h"
#include "AliEMCAL.h"
#include "TEvePointSet.h"

class Riostream;
class map;
class TTree;
class AliRun;
class AliRunLoader;
class AliEMCALLoader;
class AliEMCALGeometry;
class AliEveEMCALSModuleData;
class TClonesArray; 
class TGedFrame; 
class TEveBoxSet; 
class TEveUtil; 

class AliEveEMCALData : public TObject, public TEveRefCnt
{
 public:
  AliEveEMCALData();
  AliEveEMCALData(AliRunLoader* rl, TGeoNode* node, TGeoHMatrix* m);
  ~AliEveEMCALData();

  void SetTree(TTree* const tree);
  void SetESD(AliESDEvent* const esd);
  void SetNode(TGeoNode* const node);
  void InitEMCALGeom(AliRunLoader* const rl);
  void GetGeomInfo(Int_t id, Int_t &iSupMod, Double_t& x, Double_t& y, Double_t& z);

  void CreateAllSModules();
  void CreateSModule(Int_t sm);
  void DropAllSModules();
  void DeleteSuperModules();

  void LoadHits(TTree* const t);
  void LoadDigits(TTree* t);
  void LoadRecPoints(TTree* const t);
  void LoadHitsFromEMCALLoader(AliEMCALLoader* const emcl);
  void LoadDigitsFromEMCALLoader(AliEMCALLoader* const emcl);
  void LoadRecPointsFromEMCALLoader(AliEMCALLoader* const emcl);
  void LoadDigitsFromESD();
  void LoadRecPointsFromESD();
  void LoadRaw() const;

  AliEveEMCALSModuleData* GetSModuleData(Int_t sm);
  TEvePointSet*           GetPointSetData() const {return fPoint;};

 protected:
  AliEMCAL*         fEmcal;     // EMCal data member
  AliEMCALGeometry* fGeom;      // Data member to set/call EMCAL geometry
  TGeoNode*         fNode;      // Node for bbox definition
  TGeoHMatrix*      fHMatrix;   // matrix for local to global transformation
  TTree*            fTree;      // Tree
  AliESDEvent*      fESD;       // Esd
  Int_t             fNsm;       // Total number of Super Modules
  Int_t             fNsmfull;   // Number of full size Super Modules
  Int_t             fNsmhalf;   // Number of half size Super Modules
  std::vector<AliEveEMCALSModuleData*>   fSM;       // vector of fNsm SModules
  std::vector<AliEveEMCALSModuleData*>   fSMfull;   // vector of fNsmfull SModules
  std::vector<AliEveEMCALSModuleData*>   fSMhalf;   // vector of fNhalf SModules
  AliRunLoader*     fRunLoader; // Run Loader
  Int_t             fDebug;     // Debug option
  TEvePointSet*     fPoint;     // TEvePointSet for hits 

 private:
  AliEveEMCALData(const AliEveEMCALData &edata);            
  AliEveEMCALData& operator=(const AliEveEMCALData &edata); // Not implemented

  ClassDef(AliEveEMCALData, 0); // Base class for TRD hits visualisation
};

#endif
