// $Id$

Reve::PointSet*
tpc_hits(const char *varexp    = "TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ",
	 const char *selection = "TPC2.fArray.fR>80")
{
  // Extracts 'major' TPC hits (not the compressed ones).
  // This gives ~2.5% of all hits.

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("TPC");

  TTree* ht = rl->GetTreeH("TPC", false);

  Reve::PointSet* points = new Reve::PointSet(Form("TPC Hits '%s'", selection));
  points->SetSourceCS(TPointSelectorConsumer::TVT_RPhiZ);

  TPointSelector ps(ht, points, varexp, selection);
  ps.Select();

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerSize(2);
  points->SetMarkerColor((Color_t)3);

  gReve->AddRenderElement(points);
  gReve->Redraw3D();

  return points;
}
