// $Id$

void geom_simple()
{
  using namespace std;

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/simple_geo.root");

  Reve::GeoTopNodeRnrEl* topn_re = new Reve::GeoTopNodeRnrEl
    (gGeoManager, gGeoManager->GetTopNode());
  gReve->AddGlobalRenderElement(topn_re);
  gReve->Redraw3D();
}
