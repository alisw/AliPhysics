// $Id$

void geom_t0()
{
  gGeoManager = gEve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");

  TEveElementList* list = new TEveElementList("T0");
  gEve->AddGlobalElement(list);

  TGeoNode* node;
  TEveGeoTopNode* re;

  node = gGeoManager->GetTopVolume()->FindNode("0STR_1");
  re = new TEveGeoTopNode(gGeoManager, node);
  re->UseNodeTrans();
  gEve->AddGlobalElement(re, list);

  node = gGeoManager->GetTopVolume()->FindNode("0STL_1");
  re = new TEveGeoTopNode(gGeoManager, node);
  re->UseNodeTrans();
  gEve->AddGlobalElement(re, list);

  gEve->Redraw3D();
}
