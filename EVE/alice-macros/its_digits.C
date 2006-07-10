// $Id$

void its_digits(Int_t mode=7)
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadDigits("ITS");
  TTree* dt = rl->GetTreeD("ITS", false);

  Alieve::ITSDigitsInfo* di = new Alieve::ITSDigitsInfo();
  di->SetTree(dt);
  di->Dump();
  AliITSgeom* g = di->fGeom;

  gStyle->SetPalette(1, 0);

  gReve->DisableRedraw();

  if (mode & 1) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SPD");
    l->SetTitle("Silicon Pixel Detectors");
    l->SetMainColor((Color_t)2);
    gReve->AddRenderElement(l);
    for(Int_t i=g->GetStartSPD(); i<=g->GetLastSPD(); i++) {
      Alieve::ITSModule* m = new Alieve::ITSModule(i, di, (Color_t)2);
      gReve->AddRenderElement(l, m);
    }
  }

  if (mode & 2) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SDD");
    l->SetTitle("Silicon Drift Detectors");
    l->SetMainColor((Color_t)3);
    gReve->AddRenderElement(l);
    for(Int_t i=g->GetStartSDD(); i<=g->GetLastSDD(); i++) {
      Alieve::ITSModule* m = new Alieve::ITSModule(i, di, (Color_t)3);
      gReve->AddRenderElement(l, m);
    }
  }

  if (mode & 4) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SSD");
    l->SetTitle("Silicon Strip Detectors");
    l->SetMainColor((Color_t)4);
    gReve->AddRenderElement(l);
    for(Int_t i=g->GetStartSSD(); i<=g->GetLastSSD(); i++) {
      Alieve::ITSModule* m = new Alieve::ITSModule(i, di, (Color_t)4);
      gReve->AddRenderElement(l, m);
    } 
  }

  gReve->EnableRedraw();
}
