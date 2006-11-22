// $Id$

void geom_fmd()
{
  using namespace std;

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");

  Reve::RenderElementList* list = new Reve::RenderElementList("FMD");
  gReve->AddGlobalRenderElement(list);

  for(Int_t i=1; i<=3; ++i) {
    TGeoNode* node;
    Reve::GeoTopNodeRnrEl* re;

    //PH The line below is replaced waiting for a fix in Root
    //PH which permits to use variable siza arguments in CINT
    //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
    //PH node = gGeoManager->GetTopVolume()->FindNode(Form("F%dMT_%d", i, i));
    char form[1000];
    sprintf(form,"F%dMT_%d", i, i);
    node = gGeoManager->GetTopVolume()->FindNode(form);
    re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
    re->SetUseNodeTrans(kTRUE);
    gReve->AddGlobalRenderElement(list, re);

    sprintf(form,"F%dMB_%d", i, i);
    node = gGeoManager->GetTopVolume()->FindNode(form);
    re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
    re->SetUseNodeTrans(kTRUE);
    gReve->AddGlobalRenderElement(list, re);
  }

  gReve->Redraw3D();
}
