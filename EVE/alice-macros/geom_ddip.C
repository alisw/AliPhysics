// $Id$

void geom_ddip()
{
  using namespace std;

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("DDIP_1");

  Reve::GeoTopNodeRnrEl* its_re = 
    new Reve::GeoTopNodeRnrEl(gGeoManager, node);
  gReve->AddGlobalRenderElement(its_re);
  gReve->DrawRenderElement(its_re);
}
