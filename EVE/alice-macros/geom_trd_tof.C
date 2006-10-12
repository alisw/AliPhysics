// $Id$

void geom_trd_tof()
{
  using namespace std;

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("B077_1");

  Reve::GeoTopNodeRnrEl* its_re = 
    new Reve::GeoTopNodeRnrEl(gGeoManager, node);
  gReve->AddGlobalRenderElement(its_re);
  gReve->Redraw3D();
}
