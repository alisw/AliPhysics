// $Id$

Reve::PointSet*
tpc_hits(const char *varexp    = "TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ",
	 const char *selection = "TPC2.fArray.fR>80",
	 Option_t *option      = "goff")
{
  // Extracts 'major' TPC hits (not the compressed ones).
  // This gives ~2.5% of all hits.

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("TPC");

  TTree* ht = rl->GetTreeH("TPC", false);
  ht->SetEstimate(400*ht->GetEntries());
  ht->Draw(varexp, selection, option);
  
  Reve::PointSet* points =
    new Reve::PointSet(Form("TPC Hits '%s'", selection), ht,
		       Reve::PointSet::TVT_RPhiZ);
  points->SetTitle(Form("N=%d", points->GetN()));
  points->SetMarkerColor((Color_t)3);

  gReve->AddRenderElement(points);
  gReve->DrawRenderElement(points);

  return points;
}
