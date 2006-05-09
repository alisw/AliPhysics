// $Id$

Reve::PointSet*
its_hits(const char *varexp    = "fX:fY:fZ",
	 const char *selection = "",
	 Option_t *option      = "goff")
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("ITS");

  TTree* ht = rl->GetTreeH("ITS", false);
  ht->Draw(varexp, selection, option);
  
  Reve::PointSet* points =
    new Reve::PointSet(Form("ITS Hits '%s'", selection), ht);
  points->SetTitle(Form("N=%d", points->GetN()));
  points->SetMarkerColor((Color_t)2);

  gReve->AddRenderElement(points);
  gReve->DrawRenderElement(points);

  return points;
}
