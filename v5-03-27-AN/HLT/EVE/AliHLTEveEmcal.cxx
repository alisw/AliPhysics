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
#include "TVector3.h"
#include "AliEveHLTEventManager.h"
#include "TEveManager.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "TEveTrans.h"
#include "TGeoNode.h"
//#include "AliEMCALGeoUtils.h"
#include "AliEMCALGeometry.h"
#include "TGeoManager.h"

#include "TEveGeoShapeExtract.h"
#include "TEveGeoNode.h"

ClassImp(AliHLTEveEmcal)

AliHLTEveEmcal::AliHLTEveEmcal() : 
AliHLTEveCalo(12, "EMCAL"),
fGeoUtils(NULL)
{
  //Constructor
//  fGeoUtils = new AliEMCALGeoUtils("EMCAL_COMPLETE","EMCAL");
  fGeoUtils = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
  //fGeoUtils = new AliEMCALGeometry("EMCAL_COMPLETEV1","EMCAL");
  CreateElementList();

}


AliHLTEveEmcal::~AliHLTEveEmcal()
{
  //Destructor, not implemented
}

void AliHLTEveEmcal::CreateElementList() {
  
  TGeoNode * gEMCALNode = gGeoManager->GetTopVolume()->FindNode("XEN1_1");
  
  //gStyle->SetPalette(1, 0);
  TEveRGBAPalette* pal = new TEveRGBAPalette(0, 512);
  pal->SetLimits(0, 15);
  
  fBoxSetDigits = new TEveBoxSet[fNModules];
  fBoxSetClusters = new TEveBoxSet[fNModules];
  
  for (Int_t sm=0; sm<fNModules; ++sm) {
    
    //Digits box set
    fBoxSetDigits[sm].SetTitle(Form("Digits Module %d", sm));
    fBoxSetDigits[sm].SetName(Form("Digits Module %d", sm));
    fBoxSetDigits[sm].SetOwnIds(kTRUE);
    
    fBoxSetDigits[sm].Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
    fBoxSetDigits[sm].RefMainTrans().SetFrom(*gEMCALNode->GetDaughter(sm)->GetMatrix());
    fBoxSetDigits[sm].SetPalette(pal);
    
    AddElement(&fBoxSetDigits[sm]);
 
    //Clusters
    fBoxSetClusters[sm].SetTitle(Form("Clusters Module %d", sm));
    fBoxSetClusters[sm].SetName(Form("Clusters Module %d", sm));
    fBoxSetClusters[sm].SetOwnIds(kTRUE);
    
    fBoxSetClusters[sm].Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
    fBoxSetClusters[sm].RefMainTrans().SetFrom(*gEMCALNode->GetDaughter(sm)->GetMatrix());
    fBoxSetClusters[sm].SetPalette(pal);
    
    AddElement(&fBoxSetClusters[sm]);


  }
}

void AliHLTEveEmcal::AddClusters(Float_t * pos, Int_t module, Float_t energy) {
  //See header file for documentation

  TVector3 vec(pos);
  Int_t absId = -1;
  fGeoUtils->GetAbsCellIdFromEtaPhi(vec.Eta(), vec.Phi(), absId);
  
  TVector3 localVec;
  fGeoUtils->RelPosCellInSModule(absId, localVec);
  fBoxSetClusters[module].AddBox(15, localVec[0],  localVec[2], energy, 6.0, 6.0);
  fBoxSetClusters[module].DigitValue(static_cast<Int_t>(energy));

  //cout << "Cluster " << pos[0] << " " << pos[1] << " " << pos[2] << endl;


}

void AliHLTEveEmcal::AddDigits(UShort_t fX, UShort_t fZ, Int_t module, Float_t energy) {
  //See header file for documentation
  //Float_t x = (fX - 24)* 6.0;		
  //Float_t z = (fZ - 48) * 6.0;
	
  Int_t absid = fGeoUtils->GetAbsCellIdFromCellIndexes(module, fX, fZ);
  Double_t posX, posY, posZ;
  if(fGeoUtils->RelPosCellInSModule(absid, posX, posY, posZ)) {
    
    //cout << "digits " << posX << "  " << posY << "  " << posZ << endl;
    fBoxSetDigits[module].AddBox(15, posY, posZ, energy*10, 6.0, 6.0);

    fBoxSetDigits[module].DigitValue(static_cast<Int_t>(energy));
  } else  {
    cout <<"AliHLTEveEmcal::AddClusters:  fail"<<endl;
  }
}
