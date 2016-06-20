// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TClonesArray.h>
#include <TTree.h>
#include <TStyle.h>
#include <TGeoManager.h>
#include <TEveManager.h>
#include <TEveElement.h>

#include <AliRunLoader.h>
#include <AliTOFGeometry.h>
#include <AliEveEventManager.h>
#include <AliEveTOFDigitsInfo.h>
#include <AliEveTOFStrip.h>
#endif

void tof_digits_strips(Int_t selectedSector=-1)
{
  TClonesArray *array = 0x0;

 // Int_t nDigitsInVolume[3] = {-1, -1, -1};
  Int_t nStrips=19;
  TGeoManager *localGeoManager = AliEveEventManager::Instance()->AssertGeometry();//gEve->GetGeometry("./geometry.root");

  if (!localGeoManager) {
    printf("ERROR: no TGeo\n");
  }

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadDigits("TOF");
  TTree* dt = rl->GetTreeD("TOF", false);

  AliEveTOFDigitsInfo* di = new AliEveTOFDigitsInfo();
  di->SetTree(dt);
  di->LoadDigits();
  di->Dump();

  AliTOFGeometry* g = di->GetTOFgeometry();

  gStyle->SetPalette(1, 0);
  gEve->DisableRedraw();

  TString sPlate;
  TString bsPlate="Plate";

  Char_t sectorName[100];
  Char_t sectorTitle[200];

  TEveElementList* ll = new TEveElementList("TOF");
  ll->SetTitle("TOF detector");
  ll->SetMainColor(2);
  gEve->AddElement(ll);

  for(Int_t iSector=0; iSector<g->NSectors(); iSector++) {
    if (selectedSector!=-1 && iSector!=selectedSector) continue;
    sprintf(sectorName,"Sector%2i",iSector);
    TEveElementList* l = new TEveElementList(sectorName);
    l->SetTitle(sectorTitle);
    l->SetMainColor(2);
    gEve->AddElement(l, ll);


    for(Int_t iPlate=0; iPlate<g->NPlates(); iPlate++) {
      if ((iSector==13 || iSector==14 || iSector==15) &&
	  (iPlate==2) )continue; // staging TOF geometry with holes

      if(iPlate==2) nStrips=15;
      else nStrips=19;

      sPlate=bsPlate;
      sPlate+=iPlate;
      TEveElementList* relPlate = new TEveElementList(sPlate.Data());
      relPlate->SetMainColor(2);
      gEve->AddElement(relPlate, l);


      for(Int_t iStrip=0; iStrip<nStrips; iStrip++) {

	array = di->GetDigits(iSector,iPlate, iStrip);

	AliEveTOFStrip* m = new AliEveTOFStrip(localGeoManager,
					       iSector,iPlate,iStrip,array);
	gEve->AddElement(m, relPlate);

      }
    }
  }

  gEve->EnableRedraw();


}
