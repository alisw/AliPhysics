void tof_digits()
{  
  TClonesArray *array = 0x0;

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
  di->LoadDigits();
  di->Dump();

  AliTOFGeometry* g = di->fGeom;
 
  gStyle->SetPalette(1, 0);
  gReve->DisableRedraw();

  Reve::RenderElementList* ll = new Reve::RenderElementList("TOF");
  ll->SetTitle("TOF detector");
  ll->SetMainColor((Color_t)2);
  gReve->AddRenderElement(ll);

  for(Int_t iSector=0; iSector<g->NSectors(); iSector++) {
    
    array = di->GetDigits(iSector);
   
    Alieve::TOFSector* m = new Alieve::TOFSector(localGeoManager,iSector,array);

    gReve->AddRenderElement(m, ll);

  }

  gReve->EnableRedraw();
}
