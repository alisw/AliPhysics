// $Id$

Reve::PointSet*
fmd_hits(const char *varexp    = "fX:fY:fZ",
	 const char *selection = "")
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("FMD");

  TTree* ht = rl->GetTreeH("FMD", false);
  
  Reve::PointSet* points = new Reve::PointSet(Form("FMD Hits '%s'", selection));

  TPointSelector ps(ht, points, varexp, selection);
  ps.Select();

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerSize(.5);
  points->SetMarkerColor((Color_t)2);

  gReve->AddRenderElement(points);
  gReve->Redraw3D();

  return points;
}
