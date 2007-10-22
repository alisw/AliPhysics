// $Id$

Reve::PointSet*
its_hits(const char *varexp    = "fX:fY:fZ",
	 const char *selection = "",
         Reve::RenderElement* cont = 0)
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("ITS");

  TTree* ht = rl->GetTreeH("ITS", false);
  
  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  Reve::PointSet* points = new Reve::PointSet(Form("ITS Hits '%s'", selection));
  char form[1000];
  sprintf(form,"ITS Hits '%s'", selection);
  Reve::PointSet* points = new Reve::PointSet(form);

  TPointSelector ps(ht, points, varexp, selection);
  // ps.SetSubIdExp("fTrack:fStatus");
  ps.Select();

  if(points->Size() == 0 && gReve->GetKeepEmptyCont() == kFALSE) {
    Warning("its_hits", Form("No hits match '%s'", selection));
    delete points;
    return 0;
  }

  //PH  points->SetTitle(Form("N=%d", points->Size()));
  sprintf(form,"N=%d", points->Size());
  points->SetTitle(form);
  points->SetMarkerSize(.5);
  points->SetMarkerColor((Color_t)2);

  gReve->AddRenderElement(points, cont);
  gReve->Redraw3D();

  return points;
}
