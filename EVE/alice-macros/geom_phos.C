// $Id$

void geom_phos()
{
  using namespace std;

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");

  for(Int_t i=1; i<=5; ++i) {
    //PH The line below is replaced waiting for a fix in Root
    //PH which permits to use variable siza arguments in CINT
    //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
    //PH TGeoNode* node = gGeoManager->GetTopVolume()->FindNode(Form("PHOS_%d", i));
    char form[1000];
    sprintf(form,"PHOS_%d", i);
    TGeoNode* node = gGeoManager->GetTopVolume()->FindNode(form);

    Reve::GeoTopNodeRnrEl* re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
    re->SetUseNodeTrans(kTRUE);
    gReve->AddGlobalRenderElement(re);
  }

  gReve->Redraw3D();
}
