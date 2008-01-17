// $Id$

void geom_trd_tof()
{
  using namespace std;

  gGeoManager = gEve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("B077_1");

  TEveGeoTopNode* its_re = new TEveGeoTopNode(gGeoManager, node);
  gEve->AddGlobalElement(its_re);
  gEve->Redraw3D();
}
