// $Id$

void geom_all()
{
  gGeoManager = gEve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");

  TEveGeoTopNode* topn_re = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  gEve->AddGlobalElement(topn_re);
  gEve->Redraw3D();
}
