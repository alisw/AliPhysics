void tof_digits_strips()
{
  TClonesArray *array = 0x0;

  Int_t nDigitsInVolume[3] = {-1, -1, -1};
  Int_t nStrips=19;
  TGeoManager *localGeoManager = (TGeoManager*)gEve->GetGeometry("./geometry.root");//"$REVESYS/alice-data/alice_fullgeo.root");
  if (!localGeoManager) {
    printf("ERROR: no TGeo\n");
  }

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadDigits("TOF");
  TTree* dt = rl->GetTreeD("TOF", false);

  Alieve::TOFDigitsInfo* di = new Alieve::TOFDigitsInfo();
  di->SetTree(dt);
  di->LoadDigits();
  di->Dump();

  AliTOFGeometry* g = di->fGeom;

  gStyle->SetPalette(1, 0);
  gEve->DisableRedraw();

  TString sPlate;
  TString bsPlate="Plate";
  TString sStrip;
  TString bsStrip="Strip";
  TString sPadZ;
  TString bsPadZ="PadZ";
  TString sPadX;
  TString bsPadX="PadX";
  
  Char_t sectorName[100];
  Char_t sectorTitle[200];

  TEveElementList* ll = new TEveElementList("TOF");
  ll->SetTitle("TOF detector");
  ll->SetMainColor((Color_t)2);
  gEve->AddElement(ll);

  for(Int_t iSector=0; iSector<g->NSectors(); iSector++) {

    sprintf(sectorName,"Sector%2i",iSector);
    TEveElementList* l = new TEveElementList(sectorName);
    l->SetTitle(sectorTitle);
    l->SetMainColor((Color_t)2);
    gEve->AddElement(l, ll);


    for(Int_t iPlate=0; iPlate<g->NPlates(); iPlate++) {
      if(iPlate==2) nStrips=15;
      else nStrips=19;

      sPlate=bsPlate;
      sPlate+=iPlate;
      TEveElementList* relPlate = new TEveElementList(sPlate.Data());
      relPlate->SetMainColor((Color_t)2);
      gEve->AddElement(relPlaete, l);


      for(Int_t iStrip=0; iStrip<nStrips; iStrip++) {

	array = di->GetDigits(iSector,iPlate, iStrip);

	Alieve::TOFStrip* m = new Alieve::TOFStrip(localGeoManager,iSector,iPlate,iStrip,array);
	gEve->AddElement(m, relPlate);

      }
    }
  }

  gEve->EnableRedraw();


}
