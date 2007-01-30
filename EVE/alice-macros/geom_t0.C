// $Id$

void geom_t0()
{
  using namespace std;

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");

  Reve::RenderElementList* list = new Reve::RenderElementList("T0");
  gReve->AddGlobalRenderElement(list);

  TGeoNode* node;
  Reve::GeoTopNodeRnrEl* re;

  node = gGeoManager->GetTopVolume()->FindNode("0STR_1");
  re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
  re->SetUseNodeTrans(kTRUE);
  gReve->AddGlobalRenderElement(list, re);

  node = gGeoManager->GetTopVolume()->FindNode("0STL_1");
  re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
  re->SetUseNodeTrans(kTRUE);
  gReve->AddGlobalRenderElement(list, re);

  gReve->Redraw3D();
}
