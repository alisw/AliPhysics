// $Id$

// Load ITS digits.
// Argument mode is a bitwise or determining which layers to import:
//    1,  2 : SPD
//    4,  8 : SDD
//   16, 32 : SSD
// By default import all layers.

void its_digits(Int_t mode=63)
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

  TString sSector;
  TString bsSector="Sector";
  TString sStave;
  TString bsStave="Stave";
  TString sLadder;
  TString bsLadder="Ladder";

  Int_t i=0;
  Int_t nsec, nstave, nlad, nMod;

  if (mode & 1) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SPD0");
    l->SetTitle("SPDs' first layer");
    l->SetMainColor((Color_t)2);
    gReve->AddRenderElement(l);
    for(nsec=0; nsec<10; nsec++) {
      sSector=bsSector;
      sSector+=nsec;
      Reve::RenderElementList* relSector = new Reve::RenderElementList(sSector.Data());
      relSector->SetMainColor((Color_t)2);
      gReve->AddRenderElement(l, relSector);
      for(nstave=0; nstave<2; nstave++){
	sStave=bsStave;
	sStave += nstave;
	Reve::RenderElementList* relStave = new Reve::RenderElementList(sStave.Data());
	relStave->SetMainColor((Color_t)2);
	gReve->AddRenderElement(relSector, relStave);
	for(nMod=0; nMod<4; nMod++) {
	  Alieve::ITSModule* m = new Alieve::ITSModule(i++, di, (Color_t)2);
	  gReve->AddRenderElement(relStave, m);
	}
      }
    }
  }

  if (mode & 2) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SPD1");
    l->SetTitle("SPDs' second layer");
    l->SetMainColor((Color_t)2);
    gReve->AddRenderElement(l);

    for(nsec=0; nsec<10; nsec++) {
      sSector=bsSector;
      sSector+=nsec;
      Reve::RenderElementList* relSector = new Reve::RenderElementList(sSector.Data());
      relSector->SetMainColor((Color_t)2);
      gReve->AddRenderElement(l, relSector);
      for(nstave=0; nstave<4; nstave++){
	sStave=bsStave;
	sStave += nstave;
	Reve::RenderElementList* relStave = new Reve::RenderElementList(sStave.Data());
	relStave->SetMainColor((Color_t)2);
	gReve->AddRenderElement(relSector, relStave);
	for(nMod=0; nMod<4; nMod++) {
	  Alieve::ITSModule* m = new Alieve::ITSModule(i++, di, (Color_t)2);
	  gReve->AddRenderElement(relStave, m);
	}
      }
    }
  }

  if (mode & 4) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SDD2");
    l->SetTitle("SDDs' first layer");
    l->SetMainColor((Color_t)3);
    gReve->AddRenderElement(l);

    for(nlad=0; nlad<14; nlad++) {
      sLadder=bsLadder;
      sLadder+=nlad;
      Reve::RenderElementList* relLadder = new Reve::RenderElementList(sLadder.Data());
      relLadder->SetMainColor((Color_t)3);
      gReve->AddRenderElement(l, relLadder);
      for(nMod=0; nMod<6; nMod++) {
	Alieve::ITSModule* m = new Alieve::ITSModule(i++, di, (Color_t)3);
	gReve->AddRenderElement(relLadder, m);
      }
    }
  }

  if (mode & 8) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SDD3");
    l->SetTitle("SDDs' second layer");
    l->SetMainColor((Color_t)3);
    gReve->AddRenderElement(l);
    for(nlad=0; nlad<22; nlad++) {
      sLadder=bsLadder;
      sLadder+=nlad;
      Reve::RenderElementList* relLadder = new Reve::RenderElementList(sLadder.Data());
      relLadder->SetMainColor((Color_t)3);
      gReve->AddRenderElement(l, relLadder);
      for(nMod=0; nMod<8; nMod++) {
	Alieve::ITSModule* m = new Alieve::ITSModule(i++, di, (Color_t)3);
	gReve->AddRenderElement(relLadder, m);
      }
    }
  }

  if (mode & 16) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SSD4");
    l->SetTitle("SSDs' first layer");
    l->SetMainColor((Color_t)4);
    gReve->AddRenderElement(l);
    for(nlad=0; nlad<34; nlad++) {
      sLadder=bsLadder;
      sLadder+=nlad;
      Reve::RenderElementList* relLadder = new Reve::RenderElementList(sLadder.Data());
      relLadder->SetMainColor((Color_t)4);
      gReve->AddRenderElement(l, relLadder);
      for(nMod=0; nMod<22; nMod++) {
	Alieve::ITSModule* m = new Alieve::ITSModule(i++, di, (Color_t)4);
	gReve->AddRenderElement(relLadder, m);
      }
    }
  }

  if (mode & 32) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SSD5");
    l->SetTitle("SSDs' second layer");
    l->SetMainColor((Color_t)4);
    gReve->AddRenderElement(l);
    for(nlad=0; nlad<38; nlad++) {
      sLadder=bsLadder;
      sLadder+=nlad;
      Reve::RenderElementList* relLadder = new Reve::RenderElementList(sLadder.Data());
      relLadder->SetMainColor((Color_t)4);
      gReve->AddRenderElement(l, relLadder);
      for(nMod=0; nMod<25; nMod++) {
	Alieve::ITSModule* m = new Alieve::ITSModule(i++, di, (Color_t)4);
	gReve->AddRenderElement(relLadder, m);
      }
    }
  }

  gReve->EnableRedraw();
}
