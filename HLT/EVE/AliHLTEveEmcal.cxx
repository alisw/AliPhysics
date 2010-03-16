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
/// @brief  EMCAL class for the HLT EVE display
// Author:  Svein Lindal <slindal@fys.uio.no>

#include "AliHLTEveEmcal.h"
#include "AliHLTHOMERBlockDesc.h"
#include "TCanvas.h"
#include "AliHLTEveBase.h"
#include "TEveBoxSet.h"
#include "AliPHOSGeometry.h"
#include "TVector3.h"
#include "AliEveHOMERManager.h"
#include "TEveManager.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "TEveTrans.h"
#include "TGeoNode.h"
#include "AliEMCALGeoUtils.h"

ClassImp(AliHLTEveEmcal)

AliHLTEveEmcal::AliHLTEveEmcal() : 
AliHLTEveCalo(12, "EMCAL"),
fGeoUtils(NULL)
{
  //Constructor
  fGeoUtils = new AliEMCALGeoUtils("EMCAL_COMPLETE","EMCAL");

}


AliHLTEveEmcal::~AliHLTEveEmcal()
{
  //Destructor, not implemented
}

TEveElementList * AliHLTEveEmcal::CreateElementList() {
  
  TGeoNode * gEMCALNode = fEventManager->GetGeoManager()->GetTopVolume()->FindNode("XEN1_1");

  fElementList = new TEveElementList("EMCAL Cells");
  fElementList->SetTitle("Tooltip");
  
  //gStyle->SetPalette(1, 0);
  TEveRGBAPalette* pal = new TEveRGBAPalette(0, 512);
  pal->SetLimits(0, 15);
  
  fBoxSet = new TEveBoxSet[fNModules];
  
  for (Int_t sm=0; sm<fNModules; ++sm) {
    
    fBoxSet[sm].SetTitle(Form("Clusters Module %d", sm));
    fBoxSet[sm].SetName(Form("Clusters Module %d", sm));
    fBoxSet[sm].SetOwnIds(kTRUE);
    
    fBoxSet[sm].Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
    fBoxSet[sm].RefMainTrans().SetFrom(*gEMCALNode->GetDaughter(sm)->GetMatrix());
    fBoxSet[sm].SetPalette(pal);
    
    fElementList->AddElement(&fBoxSet[sm]);

  }
  
  return fElementList;
}

void AliHLTEveEmcal::AddClusters(Float_t * pos, Int_t module, Float_t energy) {
  //See header file for documentation

  fBoxSet[module].AddBox(15, pos[1], pos[2], energy*20, 6.1, 6.1);
  fBoxSet[module].DigitValue(static_cast<Int_t>(energy));
}

void AliHLTEveEmcal::AddDigits(UShort_t fX, UShort_t fZ, Int_t module, Float_t energy) {
  //See header file for documentation
  Int_t absid = fGeoUtils->GetAbsCellIdFromCellIndexes(module, fX, fZ);
  Double_t posX, posY, posZ;
  if(fGeoUtils->RelPosCellInSModule(absid, posX, posY, posZ)) {
    
    cout << posX << "  " << posY << "  " << posZ << endl;
    fBoxSet[module].AddBox(15, posY, posZ, energy*5, 6.1, 6.1);
    fBoxSet[module].DigitValue(static_cast<Int_t>(energy));
  } else  {
    cout <<"AliHLTEveEmcal::AddClusters:  fail"<<endl;
  }
}
