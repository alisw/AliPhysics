void its_raw(const char *input = "rawdata.root",
	     Int_t  mode       = 63,
	     Int_t  event      = 0,
	     Bool_t accumulate = kFALSE)
{
  if (gSystem->AccessPathName(input, kReadPermission))
  {
    Error("its_raw", "file '%s' not found.", input);
    return;
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

  Alieve::ITSDigitsInfo* di = new Alieve::ITSDigitsInfo();

  if (accumulate) AliLog::SetGlobalLogLevel(AliLog::kError);
  Int_t ev = 0;
  do {
    if (ev % 100 == 0) printf("Event: %d\n", ev);
    if (rawReader->NextEvent() == kFALSE)
    {
      Error("its_raw", "Reading event %d failed (requested event %d).", ev, event);
      if (accumulate)
	break;
      else
	return;
    }
    if (accumulate) di->ReadRaw(rawReader);
  } while (++ev < event);

  if ( ! accumulate) di->ReadRaw(rawReader);

  di->Dump();

  delete rawReader;
  // Could initialize ITSModule statics (?)

  AliITSgeom* g = di->fGeom;

  gStyle->SetPalette(1, 0);
  // Initialize palettes (?)

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
      sSector  = bsSector;
      sSector += nsec;
      Reve::RenderElementList* relSector = new Reve::RenderElementList(sSector.Data());
      relSector->SetMainColor((Color_t)2);
      gReve->AddRenderElement(relSector, l);
      for(nstave=0; nstave<2; nstave++){
	sStave  = bsStave;
	sStave += nstave;
	Reve::RenderElementList* relStave = new Reve::RenderElementList(sStave.Data());
	relStave->SetMainColor((Color_t)2);
	gReve->AddRenderElement(relStave, relSector);
	for(nMod=0; nMod<4; nMod++)
	{
	  if (di->GetDigits(i, 0) && di->GetDigits(i, 0)->GetEntriesFast() > 0)
	  {
	    Alieve::ITSModule* m = new Alieve::ITSModule(i, di);
	    gReve->AddRenderElement(m, relStave);
	  }
	  ++i;
	}
      }
    }
  } else {
    i += 10*2*4;
  }

  if (mode & 2) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SPD1");
    l->SetTitle("SPDs' second layer");
    l->SetMainColor((Color_t)2);
    gReve->AddRenderElement(l);

    for(nsec=0; nsec<10; nsec++) {
      sSector  = bsSector;
      sSector += nsec;
      Reve::RenderElementList* relSector = new Reve::RenderElementList(sSector.Data());
      relSector->SetMainColor((Color_t)2);
      gReve->AddRenderElement(relSector, l);
      for(nstave=0; nstave<4; nstave++){
	sStave  = bsStave;
	sStave += nstave;
	Reve::RenderElementList* relStave = new Reve::RenderElementList(sStave.Data());
	relStave->SetMainColor((Color_t)2);
	gReve->AddRenderElement(relStave, relSector);
	for(nMod=0; nMod<4; nMod++)
	{
	  if (di->GetDigits(i, 0) && di->GetDigits(i, 0)->GetEntriesFast() > 0)
	  {
	    Alieve::ITSModule* m = new Alieve::ITSModule(i, di);
	    gReve->AddRenderElement(m, relStave);
	  }
	  ++i;
	}
      }
    }
  } else {
    i += 10*4*4;
  }

  if (mode & 4) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SDD2");
    l->SetTitle("SDDs' first layer");
    l->SetMainColor((Color_t)3);
    gReve->AddRenderElement(l);

    for(nlad=0; nlad<14; nlad++) {
      sLadder  = bsLadder;
      sLadder += nlad;
      Reve::RenderElementList* relLadder = new Reve::RenderElementList(sLadder.Data());
      relLadder->SetMainColor((Color_t)3);
      gReve->AddRenderElement(relLadder, l);
      for(nMod=0; nMod<6; nMod++) {
	Alieve::ITSModule* m = new Alieve::ITSModule(i++, di);
	gReve->AddRenderElement(m, relLadder);
      }
    }
  } else {
    i += 14*6;
  }

  if (mode & 8) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SDD3");
    l->SetTitle("SDDs' second layer");
    l->SetMainColor((Color_t)3);
    gReve->AddRenderElement(l);
    for(nlad=0; nlad<22; nlad++) {
      sLadder  = bsLadder;
      sLadder += nlad;
      Reve::RenderElementList* relLadder = new Reve::RenderElementList(sLadder.Data());
      relLadder->SetMainColor((Color_t)3);
      gReve->AddRenderElement(relLadder, l);
      for(nMod=0; nMod<8; nMod++) {
	Alieve::ITSModule* m = new Alieve::ITSModule(i++, di);
	gReve->AddRenderElement(m, relLadder);
      }
    }
  } else {
    i += 22*8;
  }

  if (mode & 16) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SSD4");
    l->SetTitle("SSDs' first layer");
    l->SetMainColor((Color_t)4);
    gReve->AddRenderElement(l);
    for(nlad=0; nlad<34; nlad++) {
      sLadder  = bsLadder;
      sLadder += nlad;
      Reve::RenderElementList* relLadder = new Reve::RenderElementList(sLadder.Data());
      relLadder->SetMainColor((Color_t)4);
      gReve->AddRenderElement(relLadder, l);
      for(nMod=0; nMod<22; nMod++) {
	Alieve::ITSModule* m = new Alieve::ITSModule(i++, di);
	gReve->AddRenderElement(m, relLadder);
      }
    }
  } else {
    i += 34*22;
  }

  if (mode & 32) {
    Reve::RenderElementList* l = new Reve::RenderElementList("SSD5");
    l->SetTitle("SSDs' second layer");
    l->SetMainColor((Color_t)4);
    gReve->AddRenderElement(l);
    for(nlad=0; nlad<38; nlad++) {
      sLadder  = bsLadder;
      sLadder += nlad;
      Reve::RenderElementList* relLadder = new Reve::RenderElementList(sLadder.Data());
      relLadder->SetMainColor((Color_t)4);
      gReve->AddRenderElement(relLadder, l);
      for(nMod=0; nMod<25; nMod++) {
	Alieve::ITSModule* m = new Alieve::ITSModule(i++, di);
	gReve->AddRenderElement(m, relLadder);
      }
    }
  } else {
    i += 38*25;
  }

  gReve->EnableRedraw();
}
