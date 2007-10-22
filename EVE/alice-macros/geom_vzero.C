// $Id$

void geom_vzero()
{
  using namespace std;

  static const Reve::Exc_t eH("geom_vzero() ");

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");

  Reve::RenderElementList* list = new Reve::RenderElementList("VZero");
  gReve->AddGlobalRenderElement(list);

  TGeoNode* node = 0;
  Reve::GeoTopNodeRnrEl* re;

  TGeoNode* mnode = gGeoManager->GetTopVolume()->FindNode("VZERO_1");
  if (!mnode) {
    Error(eH, "mother node not found.");
    return;
  }

  node = mnode->GetVolume()->FindNode("V0RI_1");
  printf("opofoih %p\n", node);
  if (!node) {
    Error(eH, "V0R not found.");
    return;
  }
  re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
  re->UseNodeTrans();
  gReve->AddGlobalRenderElement(re, list);

  node = mnode->GetVolume()->FindNode("V0LE_1");
  if (!node) {
    Error(eH, "V0L not found.");
    return;
  }
  re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
  re->UseNodeTrans();
  gReve->AddGlobalRenderElement(re, list);

  gReve->Redraw3D();
}
