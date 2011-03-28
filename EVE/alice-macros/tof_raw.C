void tof_raw(Int_t newDecoder = 2)
{
  AliRawReader *rawReader = AliEveEventManager::AssertRawReader();

  TClonesArray *array = 0x0;

  TGeoManager *localGeoManager = gEve->GetDefaultGeometry();
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
