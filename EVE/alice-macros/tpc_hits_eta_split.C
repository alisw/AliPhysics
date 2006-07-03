// $Id$

void tpc_hits_eta_split(const char *varexp    =
			"TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ"
			":TPC2.fArray.Eta()",
			const char *selection = "TPC2.fArray.fR>80")
{
  // Extracts 'major' TPC hits (not the compressed ones).
  // This gives ~2.5% of all hits.

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("TPC");

  TTree* ht = rl->GetTreeH("TPC", false);

  gReve->DisableRedraw();

  Reve::PointSetArray* l = new Reve::PointSetArray("TPC hits - Eta Slices", "");
  l->SetSourceCS(TPointSelectorConsumer::TVT_RPhiZ);
  l->SetMarkerColor((Color_t)3);
  l->SetMarkerStyle(20); // full circle
  l->SetMarkerSize(2);
  
  TGListTreeItem *ti = gReve->AddRenderElement(l);

  l->InitBins(ti, "Eta", 20, -2, 2);

  TPointSelector ps(ht, l, varexp, selection);
  ps.Select();

  l->CloseBins();

  gReve->DrawRenderElement(l);
  gReve->EnableRedraw();
}
