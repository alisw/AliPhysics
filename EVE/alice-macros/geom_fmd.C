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

    node = gGeoManager->GetTopVolume()->FindNode(Form("F%dMT_%d", i, i));
    re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
    re->SetUseNodeTrans(kTRUE);
    gReve->AddGlobalRenderElement(list, re);

    node = gGeoManager->GetTopVolume()->FindNode(Form("F%dMB_%d", i, i));
    re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
    re->SetUseNodeTrans(kTRUE);
    gReve->AddGlobalRenderElement(list, re);
  }

  gReve->Redraw3D();
}
