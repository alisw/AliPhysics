#ifndef ALIEVEEMCALDATA_H
#define ALIEVEEMCALDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///
/// Fill containers for visualisation of EMCAL data structures
/// * read and store MC Hits    - read and store digits from esds or runloader
/// * read and store clusters from esds or runloader 
///
/// \author Magali Estienne <magali.estienne@cern.ch>, SUBATECH. EMCal implementation, June 2008
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS. DCal implementation + doxygen, May 2015.
///

#include <TGeoNode.h> 
#include <TGeoMatrix.h>
#include <TEveRGBAPalette.h>
#include <TEveTrans.h>
#include <TEveQuadSet.h> 
#include <TEvePointSet.h>

#include "AliESDEvent.h"
#include "AliRun.h"

class map;
class TTree;
class AliRun;
class AliRunLoader;
class AliEMCAL;
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
  
  AliEMCAL*         fEmcal;     ///< EMCal data member.
  AliEMCALGeometry* fGeom;      ///< Data member to set/call EMCAL geometry.
  TGeoNode*         fNode;      ///< Node for bbox definition.
  TGeoHMatrix*      fHMatrix;   ///< Matrix for local to global transformation.
  TTree*            fTree;      ///< Data Tree.
  AliESDEvent*      fESD;       ///< ESD event.
  
  Int_t             fNsm;       ///< Total number of Super Modules, EMCal+DCal.
  Int_t             fNsmfull;   ///< Number of full size EMCal Super-Modules.
  Int_t             fNsmhalf;   ///< Number of half size EMCal Super-Modules.
  Int_t             fNsmfullD;  ///< Number of full size DCal Super-Modules
  Int_t             fNsmhalfD;  ///< Number of half size DCal Super-Modules

  std::vector<AliEveEMCALSModuleData*>   fSM;       ///< Vector of fNsm SModules.
  std::vector<AliEveEMCALSModuleData*>   fSMfull;   ///< Vector of fNsmfull SModules.
  std::vector<AliEveEMCALSModuleData*>   fSMhalf;   ///< Vector of fNhalf SModules.
  std::vector<AliEveEMCALSModuleData*>   fSMfullD;  ///< Vector of fNsmfullD SModules.
  std::vector<AliEveEMCALSModuleData*>   fSMhalfD;  ///< Vector of fNhalfD SModules.
  
  AliRunLoader*     fRunLoader; ///< Run Loader.
  Int_t             fDebug;     ///< Debug option.
  TEvePointSet*     fPoint;     ///< TEvePointSet for hits.

 private:
  
  /// Copy constructor not implemented.
  AliEveEMCALData           (const AliEveEMCALData &edata);  
  
  /// Assignment operator not implemented.
  AliEveEMCALData& operator=(const AliEveEMCALData &edata); // Not implemented

  /// \cond CLASSIMP
  ClassDef(AliEveEMCALData, 2) ; 
  /// \endcond
  
};

#endif  //ALIEVEEMCALDATA_H
