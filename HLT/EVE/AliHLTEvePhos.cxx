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
#include "AliEveHOMERManager.h"
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
}

AliHLTEvePhos::~AliHLTEvePhos()
{
  //Destructor
  if(fGeoUtils)
    delete fGeoUtils;
  fGeoUtils = NULL;
}


TEveElementList * AliHLTEvePhos::CreateElementList() {
  //See header file for documentation

  TEveElementList * elementList  = new TEveElementList("PHOS ");
  fBoxSet = new TEveBoxSet[fNModules];

  AliPHOSGeometry * geo = AliPHOSGeometry::GetInstance("IHEP", "IHEP");
  
  TVector3 center;
  Float_t angle;
  
  // -- Create boxsets
  for(int im = 0; im < fNModules; im++) {
    
    TEveRGBAPalette* pal = new TEveRGBAPalette(0,512);
    pal->SetLimits(-1, 6);
    fBoxSet[im].SetTitle(Form("Clusters Module %d", im));
    fBoxSet[im].SetName(Form("Clusters Module %d", im));
    fBoxSet[im].SetPalette(pal);
    fBoxSet[im].Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
    fBoxSet[im].SetOwnIds(kTRUE);
    
    
    geo->GetModuleCenter(center, "CPV", im+1);
    angle = geo->GetPHOSAngle(im+1)*TMath::Pi()/180;
    
    fBoxSet[im].RefitPlex();
    TEveTrans& t = fBoxSet[im].RefMainTrans();
    t.SetupRotation(1, 2, angle );
    t.SetPos(center.X(), center.Y(), center.Z());
    
    elementList->AddElement(&fBoxSet[im]);
  }
  return elementList;
}

void AliHLTEvePhos::AddDigits(UShort_t fX, UShort_t fZ, Int_t module, Float_t energy) {
  //See header file for documentation
  Float_t x = (fX - 32)* 2.2;
  Float_t z = (fZ - 28) * 2.2;
  fBoxSet[4-module].AddBox(x, 0, z, 2.2, energy*20, 2.2);
  fBoxSet[4-module].DigitValue(static_cast<Int_t>(energy));
}


void AliHLTEvePhos::AddClusters(Float_t * pos, Int_t module, Float_t energy) {

  

  TVector3 localVector;
  TVector3 globalVector(pos);


  fGeoUtils->Global2Local(localVector, globalVector, 5-module);
  

 //See header file for documentation
  fBoxSet[4-module].AddBox(localVector.X(), -70, localVector.Z(), 2.2, -energy*20, 2.2);
  fBoxSet[4-module].DigitValue(static_cast<Int_t>(energy));
}

