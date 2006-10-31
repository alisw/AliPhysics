// $Id$

void geom_vzero()
{
  using namespace std;

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");

  Reve::RenderElementList* list = new Reve::RenderElementList("VZero");
  gReve->AddGlobalRenderElement(list);

  TGeoNode* node;
  Reve::GeoTopNodeRnrEl* re;

  node = gGeoManager->GetTopVolume()->FindNode("V0RI_1");
  re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
  re->SetUseNodeTrans(kTRUE);
  gReve->AddGlobalRenderElement(list, re);

  node = gGeoManager->GetTopVolume()->FindNode("V0LE_1");
  re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
  re->SetUseNodeTrans(kTRUE);
  gReve->AddGlobalRenderElement(list, re);

  gReve->Redraw3D();
}
