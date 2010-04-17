//*************************************************************************
// EMCAL event display
// Store the data related to each Super Module (SM)
// Possible storage of hits, digits and clusters per SM
// For visualization
// in alieve
//
//  Author: Magali Estienne (magali.estienne@cern.ch)
//  June 30 2008
//*************************************************************************

#include <Riostream.h>

#include "AliEveEMCALSModuleData.h"

#include <AliEMCALGeometry.h>

#include <TVector2.h>
#include <TVectorT.h>
#include <TClonesArray.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>

#include <EveBase/AliEveEventManager.h>

ClassImp(AliEveEMCALSModuleData)

Float_t AliEveEMCALSModuleData::fgSModuleBigBox0 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleBigBox1 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleBigBox2 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleSmallBox0 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleSmallBox1 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleSmallBox2 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleCenter0 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleCenter1 = 0.;
Float_t AliEveEMCALSModuleData::fgSModuleCenter2 = 0.;

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
  fPhiTileSize(0), fEtaTileSize(0),
  fHitArray(0),
  fDigitArray(0),
  fClusterArray(0),
  fMatrix(0),
  fHMatrix(m)
{
  //
  // Constructor
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
  fPhiTileSize(esmdata.fPhiTileSize), fEtaTileSize(esmdata.fEtaTileSize),
  fHitArray(esmdata.fHitArray),
  fDigitArray(esmdata.fDigitArray),
  fClusterArray(esmdata.fClusterArray),
  fMatrix(esmdata.fMatrix),
  fHMatrix(esmdata.fHMatrix)
{
  //
  // Copy constructor
  //

  Init(esmdata.fNsm);

}

//______________________________________________________________________________
AliEveEMCALSModuleData::~AliEveEMCALSModuleData()
{
  //
  // Destructor
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
  // Release the SM data
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

  //
  // Initialize parameters
  //

  fNsm = 12;
  fNsmfull = 10;
  fNsmhalf = 2;

  fPhiTileSize = fGeom->GetPhiTileSize();
  fEtaTileSize = fGeom->GetPhiTileSize();

  TGeoBBox* bbbox = (TGeoBBox*) fNode->GetDaughter(0) ->GetVolume()->GetShape();
  TGeoBBox* sbbox = (TGeoBBox*) fNode->GetDaughter(10)->GetVolume()->GetShape();

  fMatrix = (TGeoMatrix*) fNode->GetDaughter(sm)->GetMatrix();

  if(sm<fNsmfull)
    {
      fgSModuleBigBox0 = bbbox->GetDX();
      fgSModuleBigBox1 = bbbox->GetDY();
      fgSModuleBigBox2 = bbbox->GetDZ();
    }
  else 
    {
      fgSModuleSmallBox0 = sbbox->GetDX();
      fgSModuleSmallBox1 = sbbox->GetDY();
      fgSModuleSmallBox2 = sbbox->GetDZ();
    }
}


// ______________________________________________________________________________
void AliEveEMCALSModuleData::RegisterDigit(Int_t AbsId, Int_t isupMod, Double_t iamp, Double_t ix, Double_t iy, Double_t iz)
{
  //
  // Add a digit to this SM
  //

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

// ______________________________________________________________________________
void AliEveEMCALSModuleData::RegisterHit(Int_t AbsId, Int_t isupMod, Double_t iamp, Double_t ix, Double_t iy, Double_t iz)
{
  //
  // Add a hit to this SM
  //

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

// ______________________________________________________________________________
void AliEveEMCALSModuleData::RegisterCluster(Int_t isupMod, Double_t iamp, Double_t ix, Double_t iy, Double_t iz)
{
  //
  // Add a cluster to this SM
  //

  std::vector<Double_t> bufClu(5);
  bufClu[0] = isupMod;
  bufClu[1] = iamp;
  bufClu[2] = ix;
  bufClu[3] = iy;
  bufClu[4] = iz;

  fClusterArray.push_back(bufClu);

  fNClusters++;

}
