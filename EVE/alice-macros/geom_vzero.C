// $Id$

void geom_vzero()
{
  using namespace std;

  static const TEveException eH("geom_vzero() ");

  gGeoManager = gEve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");

  TEveElementList* list = new TEveElementList("VZero");
  gEve->AddGlobalElement(list);

  TGeoNode* node = 0;
  TEveGeoTopNode* re;

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
  re = new TEveGeoTopNode(gGeoManager, node);
  re->UseNodeTrans();
  gEve->AddGlobalElement(re, list);

  node = mnode->GetVolume()->FindNode("V0LE_1");
  if (!node) {
    Error(eH, "V0L not found.");
    return;
  }
  re = new TEveGeoTopNode(gGeoManager, node);
  re->UseNodeTrans();
  gEve->AddGlobalElement(re, list);

  gEve->Redraw3D();
}
