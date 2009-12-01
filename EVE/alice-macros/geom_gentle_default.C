void geom_gentle_default(Bool_t register_as_global=kTRUE)
{

  TFile f("$ALICE_ROOT/EVE/alice-data/gentle_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre1 = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  if (register_as_global)
  {
    gEve->AddGlobalElement(gsre1);
  }

  TFile f("$ALICE_ROOT/EVE/alice-data/gentle_rphi_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre2 = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  TFile f("$ALICE_ROOT/EVE/alice-data/gentle_rhoz_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre3 = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  TEveElement* top = gEve->GetCurrentEvent();

  AliEveMultiView *mv = AliEveMultiView::Instance();

  mv->InitGeomGentle(gsre1, gsre2, gsre3);

  gEve->FullRedraw3D(kTRUE, kTRUE);

}
