// $Id$

TEvePointSet*
acorde_hits(const char *varexp    = "ACORDE.fX:ACORDE.fY:ACORDE.fZ",
	 const char *selection = "",
         TEveElement* cont = 0)
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("ACORDE");

  TTree* ht = rl->GetTreeH("ACORDE", false);
  
  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  TEvePointSet* points = new TEvePointSet(Form("ACORDE Hits '%s'", selection));
  char form[1000];
  sprintf(form,"ACORDE Hits '%s'", selection);
  TEvePointSet* points = new TEvePointSet(form);

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  if(points->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE)
  {
    Warning("acorde_hits", Form("No hits match '%s'", selection));
    delete points;
    return 0;
  }

  //PH  points->SetTitle(Form("N=%d", points->Size()));
  sprintf(form,"N=%d", points->Size());
  points->SetTitle(form);
  points->SetMarkerSize(.5);
  points->SetMarkerColor((Color_t)2);

  gEve->AddElement(points, cont);
  gEve->Redraw3D();

  return points;
}
