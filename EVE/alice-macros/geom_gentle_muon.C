void DrawDeep(TEveGeoShape *gsre) {
  
  for (TEveElement::List_i i = gsre->BeginChildren(); i != gsre->EndChildren(); ++i) {
    TEveGeoShape* lvl = (TEveGeoShape*) *i;
    lvl->SetRnrSelf(kFALSE);
    if (!lvl->HasChildren()) {
      lvl->SetRnrSelf(kTRUE);
      lvl->SetMainColor(3);
      lvl->SetMainTransparency(50);
    }
    DrawDeep(lvl);
  }

}

TEveGeoShape* geom_gentle_muon(Bool_t updateScene = kTRUE) {

  TFile f("$ALICE_ROOT/EVE/alice-data/gentle_geo_muon.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle MUON");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  gEve->AddGlobalElement(gsre);
  f.Close();

  gsre->SetRnrSelf(kFALSE);

  DrawDeep(gsre);

  if ( updateScene ) {
    TGLViewer* v = gEve->GetDefaultGLViewer();
    v->UpdateScene();
  }

  return gsre;

}

