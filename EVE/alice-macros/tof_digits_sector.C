void tof_digits_sector(Int_t sector=0)
{
  Int_t nDigitsInVolume[3] = {-1, -1, -1};
  Int_t nStrips=19;
  TGeoManager *localGeoManager = (TGeoManager*)gReve->GetGeometry("./geometry.root");//"$REVESYS/alice-data/alice_fullgeo.root");
  if (!localGeoManager) {
    printf("ERROR: no TGeo\n");
  }

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadDigits("TOF");
  TTree* dt = rl->GetTreeD("TOF", false);

  Alieve::TOFDigitsInfo* di = new Alieve::TOFDigitsInfo();

  di->SetTree(dt);

  AliTOFGeometry* g = di->fGeom;

  gStyle->SetPalette(1, 0);
  gReve->DisableRedraw();

  Char_t sectorName[100];
  Char_t sectorTitle[200];

  Reve::RenderElementList* ll = new Reve::RenderElementList("TOF");
  ll->SetTitle("TOF detector");
  ll->SetMainColor((Color_t)2);
  gReve->AddRenderElement(ll);

  Alieve::TOFSector* m = new Alieve::TOFSector(localGeoManager, sector, dt);
  m->SetName("Sector");
  m->SetAutoTrans(kFALSE);
  m->SetTrans();
  gReve->AddRenderElement(m, ll);

  gReve->EnableRedraw();
}
