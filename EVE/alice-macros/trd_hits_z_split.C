// $Id$

void trd_hits_z_split(const char *varexp    = "fX:fY:fZ:fZ",
		      const char *selection = "")
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("TRD");

  TTree* ht = rl->GetTreeH("TRD", false);

  Reve::PointSetArray* l = new Reve::PointSetArray("TRD hits - Z Slices", "");
  l->SetMarkerColor((Color_t)7);
  l->SetMarkerStyle(20); // full circle
  l->SetMarkerSize(.5);
  
  gReve->AddRenderElement(l);
  l->InitBins("Z", 20, -360, 360);

  TPointSelector ps(ht, l, varexp, selection);
  ps.Select();

  l->CloseBins();

  gReve->Redraw3D();
}
