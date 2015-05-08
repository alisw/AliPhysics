#ifndef ALIEVEEMCALDATA_H
#define ALIEVEEMCALDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///
/// \class AliEveEMCALData
/// \brief EMCal event display data handling
///
/// Fill containers for visualisation of EMCAL data structures
/// * read and store MC Hits information (AliEMCALHit) 
/// * read and store digits from ESDs (CaloCells) or AliRunLoader (AliEMCALDigit).
/// * read and store clusters from ESDs (CaloClusters) or AliRunLoader (AliEMCALRecPoint).
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
class TClonesArray; 
class TGedFrame; 
class TEveBoxSet; 
class TEveUtil; 
class TLorentzVector;

class AliRunLoader;

class AliEMCALGeometry;
class AliEveEMCALSModuleData;

class AliEveEMCALData : public TObject, public TEveRefCnt
{
 public:
  
  AliEveEMCALData();
  AliEveEMCALData(AliRunLoader* rl, TGeoNode* node);//, TGeoHMatrix* m);
  ~AliEveEMCALData();

  void SetESD (AliESDEvent* esd)  { fESD  = esd  ; }
  void SetNode(TGeoNode *  node)  { fNode = node ; }
  
  void InitEMCALGeom();
  void GetGeomInfo(Int_t id, Int_t &iSupMod, Double_t& x, Double_t& y, Double_t& z);

  void CreateAllSModules();
  void DeleteSuperModules();

  void LoadHits();
  void LoadDigits();
  void LoadRecPoints();
  
  void LoadHitsFromEMCALLoader();
  void LoadDigitsFromEMCALLoader();
  void LoadRecPointsFromEMCALLoader();
  
  void LoadDigitsFromESD();
  void LoadRecPointsFromESD();
  
  /// Not implemented
  void LoadRaw() const { ; } 

  AliEveEMCALSModuleData* GetSModuleData(Int_t sm);
  //  TEvePointSet*           GetPointSetData() const { return fPoint ; }

 protected:
  
  AliEMCALGeometry* fGeom;      ///< Data member to set/call EMCAL geometry.
  TGeoNode*         fNode;      ///< Node for bbox definition.
//  TGeoHMatrix*      fHMatrix;   ///< Matrix for local to global transformation.
  AliESDEvent*      fESD;       ///< ESD event.
  AliRunLoader*     fRunLoader; ///< run loader
  Int_t             fNsm;       ///< Total number of Super Modules, EMCal+DCal.

  std::vector<AliEveEMCALSModuleData*>   fSM;  ///< Vector of fNsm SModules.

  //  TEvePointSet*     fPoint;     ///< TEvePointSet for hits.

  // Temporary data members
  TLorentzVector    fClusterMom ; ///< Cluster momentum

 private:
  
  AliEveEMCALData           (const AliEveEMCALData &edata);  
  
  /// Assignment operator not implemented.
  AliEveEMCALData& operator=(const AliEveEMCALData &edata); 

  /// \cond CLASSIMP
  ClassDef(AliEveEMCALData, 2) ; 
  /// \endcond
  
};

#endif  //ALIEVEEMCALDATA_H
