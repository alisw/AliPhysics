/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TClonesArray.h>
#include <TStyle.h>
#include <TGeoManager.h>
#include <TEveManager.h>
#include <TEveElement.h>

#include <AliTOFGeometry.h>
#include <AliEveEventManager.h>
#include <AliEveTOFDigitsInfo.h>
#include <AliEveTOFSector.h>
#endif

void tof_raw(Int_t newDecoder = 2)
{
    printf("*** RAW TOF ***");
    
  AliRawReader *rawReader = AliEveEventManager::AssertRawReader();

  TClonesArray *array = 0x0;

  TGeoManager *localGeoManager = AliEveEventManager::Instance()->AssertGeometry();
  if (!localGeoManager) {
    printf("ERROR: no TGeo\n");
  }

  AliEveTOFDigitsInfo* di = new AliEveTOFDigitsInfo();
  di->ReadRaw(rawReader, newDecoder);

  AliTOFGeometry* g = new AliTOFGeometry();
 
  gStyle->SetPalette(1, 0);
  gEve->DisableRedraw();

  TEveElementList* ll = new TEveElementList("TOF");
  ll->SetTitle("TOF detector");
  ll->SetMainColor(2);
  gEve->AddElement(ll);

  for(Int_t iSector=0; iSector<g->NSectors(); iSector++) {

    array = di->GetDigits(iSector);

    AliEveTOFSector* m = new AliEveTOFSector(localGeoManager,iSector,array);

    gEve->AddElement(m, ll);

  }

  delete di;
  delete g;

  gEve->EnableRedraw();
}
