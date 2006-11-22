// $Id$

Reve::PointSet*
trd_hits(const char *varexp    = "fX:fY:fZ",
	 const char *selection = "")
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("TRD");

  TTree* ht = rl->GetTreeH("TRD", false);
  
  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  Reve::PointSet* points = new Reve::PointSet(Form("TRD Hits '%s'", selection));
  char form[1000];
  sprintf(form,"TRD Hits '%s'", selection);
  Reve::PointSet* points = new Reve::PointSet(form);

  TPointSelector ps(ht, points, varexp, selection);
  ps.Select();

  //PH  points->SetTitle(Form("N=%d", points->Size()));
  sprintf(form,"N=%d", points->Size());
  points->SetTitle(form);
  points->SetMarkerSize(.5);
  points->SetMarkerColor((Color_t)7);

  gReve->AddRenderElement(points);
  gReve->Redraw3D();

  return points;
}
