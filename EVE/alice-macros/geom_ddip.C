// $Id$

void geom_ddip()
{
  using namespace std;

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("DDIP_1");

  Reve::GeoTopNodeRnrEl* re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
  re->SetUseNodeTrans(kTRUE);
  gReve->AddGlobalRenderElement(re);
  gReve->Redraw3D();
}
