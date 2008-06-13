void tof_raw(const char *input = "raw.root",
	     const char *geometry = "geometry.root",
	     Int_t  event = 0,
	     Bool_t newDecoder = kFALSE)
{

  TClonesArray *array = 0x0;

  if (gSystem->AccessPathName(input, kReadPermission))
  {
    Error("tof_raw", "file '%s' not found.", input);
    return;
  }

  TGeoManager *localGeoManager = (TGeoManager*)gEve->GetGeometry(geometry);
  if (!localGeoManager) {
    printf("ERROR: no TGeo\n");
  }

  AliRawReader *rawReader = NULL;
  TString fileName(input);
  if (fileName.EndsWith("/")) {
    rawReader = new AliRawReaderFile(fileName);
  } else if (fileName.EndsWith(".root")) {
    rawReader = new AliRawReaderRoot(fileName);
  } else if (!fileName.IsNull()) {
    rawReader = new AliRawReaderDate(fileName);
    rawReader->SelectEvents(7);
  }

  AliEveTOFDigitsInfo* di = new AliEveTOFDigitsInfo();

  for (Int_t ev=0; rawReader->NextEvent(); ev++) {
    if (ev==event) {

      if (di) {
	di->Delete();
	di = new AliEveTOFDigitsInfo();
      }

      di->ReadRaw(rawReader, newDecoder);
      continue;
    }

    else continue;

  }

  AliTOFGeometry* g = new AliTOFGeometry();
 
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
