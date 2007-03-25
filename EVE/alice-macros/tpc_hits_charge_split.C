// $Id$

void tpc_hits_charge_split(const char *varexp    =
			"TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ"
			":log(TPC2.fArray.fCharge)",
			const char *selection = "TPC2.fArray.fR>80")
{
  // Extracts 'major' TPC hits (not the compressed ones).
  // This gives ~2.5% of all hits.

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("TPC");

  TTree* ht = rl->GetTreeH("TPC", false);

  Reve::PointSetArray* l = new Reve::PointSetArray("TPC hits - Log-Charge Slices", "");
  l->SetSourceCS(TPointSelectorConsumer::TVT_RPhiZ);
  l->SetMarkerColor((Color_t)3);
  l->SetMarkerStyle(20); // full circle
  l->SetMarkerSize(.5);
  
  gReve->AddRenderElement(l);
  l->InitBins("Log Charge", 20, 0, 5);

  TPointSelector ps(ht, l, varexp, selection);
  ps.Select();

  l->CloseBins();

  gReve->Redraw3D();
}
