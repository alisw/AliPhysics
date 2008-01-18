// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/
void tof_digits()
{  
  TClonesArray *array = 0x0;

  Int_t nDigitsInVolume[3] = {-1, -1, -1};
  Int_t nStrips=19;
  TGeoManager *localGeoManager = (TGeoManager*)gEve->GetGeometry("./geometry.root");//"$REVESYS/alice-data/alice_fullgeo.root");
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

  AliTOFGeometry* g = di->fGeom;
 
  gStyle->SetPalette(1, 0);
  gEve->DisableRedraw();

  TEveElementList* ll = new TEveElementList("TOF");
  ll->SetTitle("TOF detector");
  ll->SetMainColor((Color_t)2);
  gEve->AddElement(ll);

  for(Int_t iSector=0; iSector<g->NSectors(); iSector++) {
    
    array = di->GetDigits(iSector);
   
    AliEveTOFSector* m = new AliEveTOFSector(localGeoManager,iSector,array);

    gEve->AddElement(m, ll);

  }

  gEve->EnableRedraw();
}
