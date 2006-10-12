// $Id$

Reve::PointSet*
trd_hits(const char *varexp    = "fX:fY:fZ",
	 const char *selection = "")
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("TRD");

  TTree* ht = rl->GetTreeH("TRD", false);
  
  Reve::PointSet* points = new Reve::PointSet(Form("TRD Hits '%s'", selection));

  TPointSelector ps(ht, points, varexp, selection);
  ps.Select();

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerSize(.5);
  points->SetMarkerColor((Color_t)7);

  gReve->AddRenderElement(points);
  gReve->Redraw3D();

  return points;
}
