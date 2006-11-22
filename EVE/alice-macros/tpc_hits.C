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

  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  Reve::PointSet* points = new Reve::PointSet(Form("TPC Hits '%s'", selection));
  char form[1000];
  sprintf(form,"TPC Hits '%s'", selection);
  Reve::PointSet* points = new Reve::PointSet(form);
  points->SetSourceCS(TPointSelectorConsumer::TVT_RPhiZ);

  TPointSelector ps(ht, points, varexp, selection);
  ps.Select();

  //PH  points->SetTitle(Form("N=%d", points->Size()));
  sprintf(form,"N=%d", points->Size());
  points->SetTitle(form);
  points->SetMarkerSize(.5);
  points->SetMarkerColor((Color_t)3);

  gReve->AddRenderElement(points);
  gReve->Redraw3D();

  return points;
}
