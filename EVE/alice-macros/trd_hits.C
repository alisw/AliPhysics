// $Id$

TEvePointSet*
trd_hits(const char *varexp    = "fX:fY:fZ",
	 const char *selection = "",
	 TEveElement* cont = 0)
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("TRD");

  TTree* ht = rl->GetTreeH("TRD", false);
  
  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  TEvePointSet* points = new TEvePointSet(Form("TRD Hits '%s'", selection));
  char form[1000];
  sprintf(form,"TRD Hits '%s'", selection);
  TEvePointSet* points = new TEvePointSet(form);

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  if (points->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("trd_hits", Form("No hits match '%s'", selection));
    delete points;
    return 0;
  }

  //PH  points->SetTitle(Form("N=%d", points->Size()));
  sprintf(form,"N=%d", points->Size());
  points->SetTitle(form);
  points->SetMarkerSize(.5);
  points->SetMarkerColor((Color_t)7);

  gEve->AddElement(points, cont);
  gEve->Redraw3D();

  return points;
}
