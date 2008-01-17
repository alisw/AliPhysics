// $Id$

void its_hits_layer_split(const char *varexp    = "fX:fY:fZ:GetLayer()",
                          const char *selection = "")
{
  // Extracts 'major' TPC hits (not the compressed ones).
  // This gives ~2.5% of all hits.

  printf("THIS SCRIPT DOES NOT WORK.\n"
	 "GetLayer() crashes when trying to load ITS geometry.\n"
	 "Needs to be fixed together with ITS experts.\n");
  return;

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("ITS");

  TTree* ht = rl->GetTreeH("ITS", false);

  TEvePointSetArray* l = new TEvePointSetArray("ITS hits - Layer Slices", "");
  l->SetMarkerColor((Color_t)2);
  l->SetMarkerStyle(2); // cross
  l->SetMarkerSize(.2);

  gEve->AddElement(l);
  l->InitBins("Layer", 6, 0.5, 6.5);

  TEvePointSelector ps(ht, l, varexp, selection);
  ps.Select();

  l->CloseBins();

  gEve->Redraw3D();
}
