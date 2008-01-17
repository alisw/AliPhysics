// $Id$

TEvePointSet*
vzero_hits(const char *varexp    = "fX:fY:fZ",
	   const char *selection = "")
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("VZERO");

  TTree* ht = rl->GetTreeH("VZERO", false);
  
  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  TEvePointSet* points = new TEvePointSet(Form("VZERO Hits '%s'", selection));
  char form[1000];
  sprintf(form,"VZERO Hits '%s'", selection);
  TEvePointSet* points = new TEvePointSet(form);

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  //PH  points->SetTitle(Form("N=%d", points->Size()));
  sprintf(form,"N=%d", points->Size());
  points->SetTitle(form);
  points->SetMarkerSize(.5);
  points->SetMarkerColor((Color_t)2);

  gEve->AddElement(points);
  gEve->Redraw3D();

  return points;
}
