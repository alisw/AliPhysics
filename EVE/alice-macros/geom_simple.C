// $Id$

void geom_simple()
{
  gGeoManager = gEve->GetGeometry("$REVESYS/alice-data/simple_geo.root");

  TEveGeoTopNode* topn_re = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  gEve->AddGlobalElement(topn_re);
  gEve->Redraw3D();
}
