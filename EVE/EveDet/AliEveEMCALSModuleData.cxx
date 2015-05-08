/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TGeoBBox.h>

#include "AliEMCALGeometry.h"

#include "AliEveEMCALSModuleData.h"

class TClonesArray;
class TGeoNode;
//class TGeoMatrix;
class TVector2;

class AliEveEventManager;

/// \cond CLASSIMP
ClassImp(AliEveEMCALSModuleData) ;
/// \endcond

Float_t AliEveEMCALSModuleData::fgSModuleBigBox0 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleBigBox1 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleBigBox2 = 0.;

Float_t AliEveEMCALSModuleData::fgSModuleSmallBox0 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleSmallBox1 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleSmallBox2 = 0.;

Float_t AliEveEMCALSModuleData::fgSModuleDCalBox0 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleDCalBox1 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleDCalBox2 = 0.;

Float_t AliEveEMCALSModuleData::fgSModuleSmallDBox0 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleSmallDBox1 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleSmallDBox2 = 0.;

//
// Constructor
//
//______________________________________________________________________________
AliEveEMCALSModuleData::AliEveEMCALSModuleData(Int_t sm,AliEMCALGeometry* geom, TGeoNode* node): //, TGeoHMatrix* m) :
  TObject(),
  fGeom(geom),
  fNode(node),
  fSmId(sm),
  fNsm(0),
  fNDigits(0),
  fNClusters(0),
  fNHits(0),
  fPhiTileSize(0), fEtaTileSize(0),
  fHitArray(0),
  fDigitArray(0),
  fClusterArray(0)
//  fMatrix(0), 
//  fHMatrix(m) 
{
  Init(sm);
}

///
/// Copy constructor
///
//______________________________________________________________________________
  AliEveEMCALSModuleData::AliEveEMCALSModuleData(const AliEveEMCALSModuleData &esmdata) :
  TObject(),
  fGeom(esmdata.fGeom),
  fNode(esmdata.fNode),
  fSmId(esmdata.fSmId),
  fNsm(esmdata.fNsm),
  fNDigits(esmdata.fNDigits),
  fNClusters(esmdata.fNClusters),
  fNHits(esmdata.fNHits),
  fPhiTileSize(esmdata.fPhiTileSize), fEtaTileSize(esmdata.fEtaTileSize),
  fHitArray(esmdata.fHitArray),
  fDigitArray(esmdata.fDigitArray),
  fClusterArray(esmdata.fClusterArray)
//  fMatrix(esmdata.fMatrix),
//  fHMatrix(esmdata.fHMatrix)
{
  Init(esmdata.fNsm);
}

///
/// Destructor
///
//______________________________________________________________________________
AliEveEMCALSModuleData::~AliEveEMCALSModuleData()
{
  if(!fHitArray.empty())
    fHitArray.clear();
  
  if(!fDigitArray.empty())
    fDigitArray.clear();
  
  if(!fClusterArray.empty())
    fClusterArray.clear();
}

///
/// Release the SM data.
///
//______________________________________________________________________________
void AliEveEMCALSModuleData::DropData()
{
  fNDigits   = 0;
  fNClusters = 0;
  fNHits     = 0;

  if(!fHitArray.empty())
    fHitArray.clear();

  if(!fDigitArray.empty())
    fDigitArray.clear();

  if(!fClusterArray.empty())
    fClusterArray.clear();

  return;
}

///
/// Initialize parameters
///
// ______________________________________________________________________________
void AliEveEMCALSModuleData::Init(Int_t sm)
{
  fNsm         = fGeom->GetNumberOfSuperModules();
  fPhiTileSize = fGeom->GetPhiTileSize();
  fEtaTileSize = fGeom->GetPhiTileSize();

  //fMatrix = (TGeoMatrix*) fNode->GetDaughter(sm)->GetMatrix();

  TGeoBBox * bbox  = (TGeoBBox*) fNode->GetDaughter(sm)->GetVolume()->GetShape();
  
  if(sm < 10)
  {
    fgSModuleBigBox0    = bbox->GetDX();
    fgSModuleBigBox1    = bbox->GetDY();
    fgSModuleBigBox2    = bbox->GetDZ();
  }
  else if(sm < 12) 
  {
    fgSModuleSmallBox0  = bbox->GetDX();
    fgSModuleSmallBox1  = bbox->GetDY();
    fgSModuleSmallBox2  = bbox->GetDZ();
  }  
  else if(sm < 18) 
  {
    fgSModuleDCalBox0   = bbox->GetDX();
    fgSModuleDCalBox1   = bbox->GetDY();
    fgSModuleDCalBox2   = bbox->GetDZ();
  }  
  else if(sm < 20) 
  {
    fgSModuleSmallDBox0 = bbox->GetDX();
    fgSModuleSmallDBox1 = bbox->GetDY();
    fgSModuleSmallDBox2 = bbox->GetDZ();
  }
}

///
/// Add a digit to this SM
///
// ______________________________________________________________________________
void AliEveEMCALSModuleData::RegisterDigit(Int_t AbsId, Int_t isupMod, Double_t iamp, Double_t ix, Double_t iy, Double_t iz)
{
  std::vector<Double_t> bufDig(6);
  bufDig[0] = AbsId;
  bufDig[1] = isupMod;
  bufDig[2] = iamp;
  bufDig[3] = ix;
  bufDig[4] = iy;
  bufDig[5] = iz;

  fDigitArray.push_back(bufDig);

  fNDigits++;
}

///
/// Add a hit to this SM
///
// ______________________________________________________________________________
void AliEveEMCALSModuleData::RegisterHit(Int_t AbsId, Int_t isupMod, Double_t iamp, Double_t ix, Double_t iy, Double_t iz)
{
  std::vector<Float_t> bufHit(6);
  bufHit[0] = AbsId;
  bufHit[1] = isupMod;
  bufHit[2] = iamp;
  bufHit[3] = ix;
  bufHit[4] = iy;
  bufHit[5] = iz;

  fHitArray.push_back(bufHit);

  fNHits++;
}

///
/// Add a cluster to this SM
///
// ______________________________________________________________________________
void AliEveEMCALSModuleData::RegisterCluster(Int_t isupMod, Double_t iamp, Double_t ix, Double_t iy, Double_t iz)
{
  std::vector<Double_t> bufClu(5);
  bufClu[0] = isupMod;
  bufClu[1] = iamp;
  bufClu[2] = ix;
  bufClu[3] = iy;
  bufClu[4] = iz;

  fClusterArray.push_back(bufClu);

  fNClusters++;
}
