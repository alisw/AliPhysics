// $Id$

void geom_tpc()
{
  using namespace std;

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("TPC_M_1");

  Reve::GeoTopNodeRnrEl* tpc_re = 
    new Reve::GeoTopNodeRnrEl(gGeoManager, node);
  gReve->AddGlobalRenderElement(tpc_re);
  gReve->DrawRenderElement(tpc_re);
}
