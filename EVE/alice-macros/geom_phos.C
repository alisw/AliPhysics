// $Id$

void geom_phos()
{
  gGeoManager = gEve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");

  TEveElementList* list = new TEveElementList("PHOS");
  gEve->AddGlobalElement(list);

  for(Int_t i=1; i<=5; ++i) {
    //PH The line below is replaced waiting for a fix in Root
    //PH which permits to use variable siza arguments in CINT
    //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
    //PH TGeoNode* node = gGeoManager->GetTopVolume()->FindNode(Form("PHOS_%d", i));
    char form[1000];
    sprintf(form,"PHOS_%d", i);
    TGeoNode* node = gGeoManager->GetTopVolume()->FindNode(form);

    TEveGeoTopNode* re = new TEveGeoTopNode(gGeoManager, node);
    re->UseNodeTrans();
    gEve->AddGlobalElement(re, list);
  }

  gEve->Redraw3D();
}
