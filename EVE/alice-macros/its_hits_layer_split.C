// $Id$

void its_hits_layer_split(const char *varexp    = "fX:fY:fZ:fLayer",
                          const char *selection = "")
{
  // Extracts 'major' TPC hits (not the compressed ones).
  // This gives ~2.5% of all hits.

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("ITS");

  TTree* ht = rl->GetTreeH("ITS", false);

  Reve::PointSetArray* l = new Reve::PointSetArray("ITS hits - Layer Slices", "");
  l->SetMarkerColor((Color_t)2);
  l->SetMarkerStyle(2); // cross
  l->SetMarkerSize(.2);

  gReve->AddRenderElement(l);
  l->InitBins("Layer", 6, 0.5, 6.5);

  TPointSelector ps(ht, l, varexp, selection);
  ps.Select();

  l->CloseBins();

  gReve->Redraw3D();
}
