// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TStyle.h>
#include <TGeoManager.h>
#include <TEveManager.h>
#include <TEveElement.h>

#include <AliRunLoader.h>
#include <AliTOFGeometry.h>
#include <AliEveEventManager.h>
#include <AliEveTOFDigitsInfo.h>
#include <AliEveTOFSector.h>
#endif

void tof_digits_sector(Int_t sector=0)
{
  
  TGeoManager *localGeoManager = AliEveEventManager::Instance()->AssertGeometry(); //gEve->GetGeometry("./geometry.root");
  if (!localGeoManager) {
    printf("ERROR: no TGeo\n");
  }

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadDigits("TOF");
  TTree* dt = rl->GetTreeD("TOF", false);

  AliEveTOFDigitsInfo* di = new AliEveTOFDigitsInfo();

  di->SetTree(dt);

  //AliTOFGeometry* g = di->GetTOFgeometry();

  gStyle->SetPalette(1, 0);
  gEve->DisableRedraw();

  //Char_t sectorName[100];
  //Char_t sectorTitle[200];

  TEveElementList* ll = new TEveElementList("TOF");
  ll->SetTitle("TOF detector");
  ll->SetMainColor(2);
  gEve->AddElement(ll);

  AliEveTOFSector* m = new AliEveTOFSector(localGeoManager, sector, dt);
  m->SetName("Sector");
  m->SetAutoTrans(kFALSE);
  m->SetTrans();
  gEve->AddElement(m, ll);

  gEve->EnableRedraw();
}
