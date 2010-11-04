/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Svein Lindal <slindal@fys.uio.no   >                  *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTEvePhos.cxx
/// @author Svein Lindal <slindal@fys.uio.no>
/// @brief  PHOS class for the HLT EVE display

#include "AliHLTEvePhos.h"
#include "AliHLTHOMERBlockDesc.h"
#include "TCanvas.h"
#include "TEveBoxSet.h"
#include "AliPHOSGeoUtils.h"
#include "AliPHOSGeometry.h"
#include "TVector3.h"
#include "AliEveHLTEventManager.h"
#include "TEveManager.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "TEveTrans.h"

ClassImp(AliHLTEvePhos)

AliHLTEvePhos::AliHLTEvePhos() : 
AliHLTEveCalo(5, "PHOS"),
fGeoUtils(NULL)
{
  // Constructor.
  fGeoUtils = new AliPHOSGeoUtils("PHOS", "noCPV");
  CreateElementList();

}

AliHLTEvePhos::~AliHLTEvePhos()
{
  //Destructor
  if(fGeoUtils)
    delete fGeoUtils;
  fGeoUtils = NULL;
}

void AliHLTEvePhos::CreateElementList() {
  //See header file for documentation

  fBoxSetClusters = new TEveBoxSet[fNModules];
  fBoxSetDigits = new TEveBoxSet[fNModules];

  AliPHOSGeometry * geo = AliPHOSGeometry::GetInstance("IHEP", "IHEP");
  
  TVector3 center;
  Float_t angle;
  
  // -- Create boxsets
  for(int im = 0; im < fNModules; im++) {
    
    TEveRGBAPalette* pal = new TEveRGBAPalette(0,512);
    pal->SetLimits(-1, 45);

    //Create clusters box set
    fBoxSetClusters[im].SetTitle(Form("Clusters Module %d", im));
    fBoxSetClusters[im].SetName(Form("Clusters Module %d", im));
    fBoxSetClusters[im].SetPalette(pal);
    fBoxSetClusters[im].Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
    fBoxSetClusters[im].SetOwnIds(kTRUE);
    
    
    geo->GetModuleCenter(center, "CPV", im+1);
    angle = geo->GetPHOSAngle(im+1)*TMath::Pi()/180;
    
    fBoxSetClusters[im].RefitPlex();
    TEveTrans& t = fBoxSetClusters[im].RefMainTrans();
    t.SetupRotation(1, 2, angle );
    t.SetPos(center.X(), center.Y(), center.Z());
    
    AddElement(&fBoxSetClusters[im]);


    //Create digits box set
    fBoxSetDigits[im].SetTitle(Form("Digits Module %d", im));
    fBoxSetDigits[im].SetName(Form("Digits Module %d", im));
    fBoxSetDigits[im].SetPalette(pal);
    fBoxSetDigits[im].Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
    fBoxSetDigits[im].SetOwnIds(kTRUE);
    
    
    geo->GetModuleCenter(center, "CPV", im+1);
    angle = geo->GetPHOSAngle(im+1)*TMath::Pi()/180;
    
    fBoxSetDigits[im].RefitPlex();
    TEveTrans& t2 = fBoxSetDigits[im].RefMainTrans();
    t2.SetupRotation(1, 2, angle );
    t2.SetPos(center.X(), center.Y(), center.Z());
    
    AddElement(&fBoxSetDigits[im]);


  }

}

void AliHLTEvePhos::AddDigits(UShort_t fX, UShort_t fZ, Int_t module, Float_t energy) {
  //See header file for documentation
  Float_t x = (fX - 32)* 2.2;
  Float_t z = (fZ - 28) * 2.2;
  fBoxSetDigits[4-module].AddBox(x, 0, z, 2.2, energy*20, 2.2);
  fBoxSetDigits[4-module].DigitValue(static_cast<Int_t>(energy));
}

void AliHLTEvePhos::ProcessESDCluster(AliESDCaloCluster * cluster) {
  
  Float_t pos[3];
  cluster->GetPosition(pos);
  TVector3 vec = TVector3(pos);
  Float_t phi = vec.Phi() *180 / TMath::Pi();
  if(phi < 0) phi += 360;

  Int_t mid = 0;
  if(phi > 320) return;
  else if(phi > 300) mid = 4;
  else if(phi > 280) mid = 3;
  else if(phi > 260) mid = 2;
  else return;
  
  Int_t nCells = cluster->GetNCells();
  AddClusters(pos, mid, cluster->E(), nCells);

}

///____________________________________________________________________________________
void AliHLTEvePhos::AddClusters(Float_t * pos, Int_t module, Float_t energy) {

  TVector3 localVector;
  TVector3 globalVector(pos);

  fGeoUtils->Global2Local(localVector, globalVector, 5-module);

 //See header file for documentation
  fBoxSetClusters[4-module].AddBox(localVector.X(), 0, localVector.Z(), 2.2, energy*40, 2.2 );
  fBoxSetClusters[4-module].DigitValue(static_cast<Int_t>(energy*10));
}


///_____________________________________________________________________________
void AliHLTEvePhos::AddClusters(Float_t * pos, Int_t module, Float_t energy, Int_t nCells) {

  TVector3 localVector;
  TVector3 globalVector(pos);

  fGeoUtils->Global2Local(localVector, globalVector, 5-module);

  Float_t cf = 2.2*TMath::Sqrt(nCells);
 //See header file for documentation
  fBoxSetClusters[4-module].AddBox(localVector.X(), 0, localVector.Z(), cf, energy*40, cf );
  fBoxSetClusters[4-module].DigitValue(static_cast<Int_t>(energy*10));
}

