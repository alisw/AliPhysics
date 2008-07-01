//*************************************************************************
// EMCAL event display
// Store the data related to each Super Module (SM)
// Possible storage of hits, digits and clusters per SM
//
//  Author: Magali Estienne (magali.estienne@cern.ch)
//  June 30 2008
//*************************************************************************

#include <Riostream.h>
#include <vector>

#include "AliEveEMCALSModuleData.h"

#include <AliEMCALGeometry.h>

#include <TVector2.h>
#include <TVectorT.h>
#include <TClonesArray.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>

#include <EveBase/AliEveEventManager.h>

///////////////////////////////////////////////////////////////////////////////
///
/// AliEveEMCALSModuleData: geometry and digits
///
///////////////////////////////////////////////////////////////////////////////


ClassImp(AliEveEMCALSModuleData)

Float_t AliEveEMCALSModuleData::fSModuleBigBox0 = 0.;
Float_t AliEveEMCALSModuleData::fSModuleBigBox1 = 0.;
Float_t AliEveEMCALSModuleData::fSModuleBigBox2 = 0.;
Float_t AliEveEMCALSModuleData::fSModuleSmallBox0 = 0.;
Float_t AliEveEMCALSModuleData::fSModuleSmallBox1 = 0.;
Float_t AliEveEMCALSModuleData::fSModuleSmallBox2 = 0.;
Float_t AliEveEMCALSModuleData::fSModuleCenter0 = 0.;
Float_t AliEveEMCALSModuleData::fSModuleCenter1 = 0.;
Float_t AliEveEMCALSModuleData::fSModuleCenter2 = 0.;

//______________________________________________________________________________
AliEveEMCALSModuleData::AliEveEMCALSModuleData(Int_t sm,AliEMCALGeometry* geom, TGeoNode* node, TGeoHMatrix* m) :
  TObject(),
  fGeom(geom),
  fNode(node),
  fSmId(sm),
  fNsm(0),
  fNsmfull(0),
  fNsmhalf(0),
  fNDigits(0),
  fNClusters(0),
  fNHits(0),
  fHitArray(0),
  fDigitArray(0),
  fClusterArray(0),
  fMatrix(0),
  fHMatrix(m)
{
  //
  // constructor
  //

  Init(sm);

}

//______________________________________________________________________________
  AliEveEMCALSModuleData::AliEveEMCALSModuleData(const AliEveEMCALSModuleData &esmdata) :
  TObject(),
  fGeom(esmdata.fGeom),
  fNode(esmdata.fNode),
  fSmId(esmdata.fSmId),
  fNsm(esmdata.fNsm),
  fNsmfull(esmdata.fNsmfull),
  fNsmhalf(esmdata.fNsmhalf),
  fNDigits(esmdata.fNDigits),
  fNClusters(esmdata.fNClusters),
  fNHits(esmdata.fNHits),
  fHitArray(esmdata.fHitArray),
  fDigitArray(esmdata.fDigitArray),
  fClusterArray(esmdata.fClusterArray),
  fMatrix(esmdata.fMatrix),
  fHMatrix(esmdata.fHMatrix)
{
  //
  // constructor
  //

  Init(esmdata.fNsm);

}

//______________________________________________________________________________
AliEveEMCALSModuleData::~AliEveEMCALSModuleData()
{
  //
  // destructor
  //

  if(!fHitArray.empty()){
    fHitArray.clear();
  }
  if(!fDigitArray.empty()){
    fDigitArray.clear();
  }
  if(!fClusterArray.empty()){
    fClusterArray.clear();
  }

}

//______________________________________________________________________________
void AliEveEMCALSModuleData::DropData()
{
  //
  // release the sm data
  //

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

// ______________________________________________________________________________
void AliEveEMCALSModuleData::Init(Int_t sm)
{
  fNsm = 12;
  fNsmfull = 10;
  fNsmhalf = 2;

  fPhiTileSize = fGeom->GetPhiTileSize();
  fEtaTileSize = fGeom->GetPhiTileSize();

  TGeoBBox* bbbox = (TGeoBBox*) fNode->GetDaughter(0) ->GetVolume()->GetShape();
  bbbox->Dump();
  TGeoBBox* sbbox = (TGeoBBox*) fNode->GetDaughter(10)->GetVolume()->GetShape();
  sbbox->Dump();

  fMatrix = (TGeoMatrix*) fNode->GetDaughter(sm)->GetMatrix();

  if(sm<fNsmfull)
    {
      fSModuleBigBox0 = bbbox->GetDX();
      fSModuleBigBox1 = bbbox->GetDY();
      fSModuleBigBox2 = bbbox->GetDZ();
    }
  else 
    {
      fSModuleSmallBox0 = sbbox->GetDX();
      fSModuleSmallBox1 = sbbox->GetDY();
      fSModuleSmallBox2 = sbbox->GetDZ();
    }
}


// ______________________________________________________________________________
void AliEveEMCALSModuleData::RegisterDigit(Int_t AbsId, Int_t isupMod, Float_t iamp, Float_t ix, Float_t iy, Float_t iz)
{
  //
  // add a digit to this sm
  //

  vector<Float_t> bufDig(6);
  bufDig[0] = AbsId;
  bufDig[1] = isupMod;
  bufDig[2] = iamp;
  bufDig[3] = ix;
  bufDig[4] = iy;
  bufDig[5] = iz;

  cout << "bufDig[0]: " <<  bufDig[0] << ", bufDig[1]: " <<  bufDig[1] << ", bufDig[2]: " <<  bufDig[2] <<
    ", bufDig[3]: " <<  bufDig[3] << ", bufDig[4]: " <<  bufDig[4] << ", bufDig[5]: " <<  bufDig[5] << endl;

  fDigitArray.push_back(bufDig);

  fNDigits++;

}

// ______________________________________________________________________________
void AliEveEMCALSModuleData::RegisterHit(Int_t AbsId, Int_t isupMod, Float_t iamp, Float_t ix, Float_t iy, Float_t iz)
{
  //
  // add a hit to this sm
  //

  vector<Float_t> bufHit(6);
  bufHit[0] = AbsId;
  bufHit[1] = isupMod;
  bufHit[2] = iamp;
  bufHit[3] = ix;
  bufHit[4] = iy;
  bufHit[5] = iz;

  fHitArray.push_back(bufHit);

  fNHits++;

}

// ______________________________________________________________________________
void AliEveEMCALSModuleData::RegisterCluster(Int_t isupMod, Float_t iamp, Float_t ix, Float_t iy, Float_t iz)
{
  //
  // add a cluster to this sm
  //

  vector<Float_t> bufClu(5);
  bufClu[0] = isupMod;
  bufClu[1] = iamp;
  bufClu[2] = ix;
  bufClu[3] = iy;
  bufClu[4] = iz;

  fClusterArray.push_back(bufClu);

  fNClusters++;

}
