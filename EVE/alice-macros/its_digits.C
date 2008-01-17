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
  // Could initialize ITSModule statics (?)

  AliITSgeom* g = di->fGeom;

  gStyle->SetPalette(1, 0);
  // Initialize palettes (?)

  gEve->DisableRedraw();

  TString sSector;
  TString bsSector="Sector";
  TString sStave;
  TString bsStave="Stave";
  TString sLadder;
  TString bsLadder="Ladder";

  Int_t i=0;
  Int_t nsec, nstave, nlad, nMod;

  if (mode & 1) {
    TEveElementList* l = new TEveElementList("SPD0");
    l->SetTitle("SPDs' first layer");
    l->SetMainColor((Color_t)2);
    gEve->AddElement(l);
    for (nsec=0; nsec<10; nsec++) {
      sSector  = bsSector;
      sSector += nsec;
      TEveElementList* relSector = new TEveElementList(sSector.Data());
      relSector->SetMainColor((Color_t)2);
      gEve->AddElement(relSector, l);
      for (nstave=0; nstave<2; nstave++){
	sStave  = bsStave;
	sStave += nstave;
	TEveElementList* relStave = new TEveElementList(sStave.Data());
	relStave->SetMainColor((Color_t)2);
	gEve->AddElement(relStave, relSector);
	for (nMod=0; nMod<4; nMod++) {
	  Alieve::ITSModule* m = new Alieve::ITSModule(i++, di);
	  gEve->AddElement(m, relStave);
	}
      }
    }
  } else {
    i += 10*2*4;
  }

  if (mode & 2) {
    TEveElementList* l = new TEveElementList("SPD1");
    l->SetTitle("SPDs' second layer");
    l->SetMainColor((Color_t)2);
    gEve->AddElement(l);

    for (nsec=0; nsec<10; nsec++) {
      sSector  = bsSector;
      sSector += nsec;
      TEveElementList* relSector = new TEveElementList(sSector.Data());
      relSector->SetMainColor((Color_t)2);
      gEve->AddElement(relSector, l);
      for (nstave=0; nstave<4; nstave++){
	sStave  = bsStave;
	sStave += nstave;
	TEveElementList* relStave = new TEveElementList(sStave.Data());
	relStave->SetMainColor((Color_t)2);
	gEve->AddElement(relStave, relSector);
	for (nMod=0; nMod<4; nMod++) {
	  Alieve::ITSModule* m = new Alieve::ITSModule(i++, di);
	  gEve->AddElement(m, relStave);
	}
      }
    }
  } else {
    i += 10*4*4;
  }

  if (mode & 4) {
    TEveElementList* l = new TEveElementList("SDD2");
    l->SetTitle("SDDs' first layer");
    l->SetMainColor((Color_t)3);
    gEve->AddElement(l);

    for (nlad=0; nlad<14; nlad++) {
      sLadder  = bsLadder;
      sLadder += nlad;
      TEveElementList* relLadder = new TEveElementList(sLadder.Data());
      relLadder->SetMainColor((Color_t)3);
      gEve->AddElement(relLadder, l);
      for (nMod=0; nMod<6; nMod++) {
	Alieve::ITSModule* m = new Alieve::ITSModule(i++, di);
	gEve->AddElement(m, relLadder);
      }
    }
  } else {
    i += 14*6;
  }

  if (mode & 8) {
    TEveElementList* l = new TEveElementList("SDD3");
    l->SetTitle("SDDs' second layer");
    l->SetMainColor((Color_t)3);
    gEve->AddElement(l);
    for (nlad=0; nlad<22; nlad++) {
      sLadder  = bsLadder;
      sLadder += nlad;
      TEveElementList* relLadder = new TEveElementList(sLadder.Data());
      relLadder->SetMainColor((Color_t)3);
      gEve->AddElement(relLadder, l);
      for (nMod=0; nMod<8; nMod++) {
	Alieve::ITSModule* m = new Alieve::ITSModule(i++, di);
	gEve->AddElement(m, relLadder);
      }
    }
  } else {
    i += 22*8;
  }

  if (mode & 16) {
    TEveElementList* l = new TEveElementList("SSD4");
    l->SetTitle("SSDs' first layer");
    l->SetMainColor((Color_t)4);
    gEve->AddElement(l);
    for (nlad=0; nlad<34; nlad++) {
      sLadder  = bsLadder;
      sLadder += nlad;
      TEveElementList* relLadder = new TEveElementList(sLadder.Data());
      relLadder->SetMainColor((Color_t)4);
      gEve->AddElement(relLadder, l);
      for (nMod=0; nMod<22; nMod++) {
	Alieve::ITSModule* m = new Alieve::ITSModule(i++, di);
	gEve->AddElement(m, relLadder);
      }
    }
  } else {
    i += 34*22;
  }

  if (mode & 32) {
    TEveElementList* l = new TEveElementList("SSD5");
    l->SetTitle("SSDs' second layer");
    l->SetMainColor((Color_t)4);
    gEve->AddElement(l);
    for (nlad=0; nlad<38; nlad++) {
      sLadder  = bsLadder;
      sLadder += nlad;
      TEveElementList* relLadder = new TEveElementList(sLadder.Data());
      relLadder->SetMainColor((Color_t)4);
      gEve->AddElement(relLadder, l);
      for (nMod=0; nMod<25; nMod++) {
	Alieve::ITSModule* m = new Alieve::ITSModule(i++, di);
	gEve->AddElement(m, relLadder);
      }
    }
  } else {
    i += 38*25;
  }

  gEve->EnableRedraw();
}
