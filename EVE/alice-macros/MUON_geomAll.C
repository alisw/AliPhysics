// $Id$

void MUON_geomAll()
{

  using namespace std;

  TString dataPath = TString(Alieve::gEvent->GetTitle());
  dataPath.Append("/geometry.root");
  gGeoManager = gReve->GetGeometry(dataPath.Data());

  Reve::GeoTopNodeRnrEl* topn_re = new Reve::GeoTopNodeRnrEl
    (gGeoManager, gGeoManager->GetTopNode());
  
  gReve->AddGlobalRenderElement(topn_re);

  gReve->Redraw3D(kTRUE);
  
}
