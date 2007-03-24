// $Id$

void MUON_geom()
{
  using namespace std;

  TString dataPath = TString(Alieve::gEvent->GetTitle());
  dataPath.Append("/geometry.root");
  gGeoManager = gReve->GetGeometry(dataPath.Data());

  TGeoNode *node1 = gGeoManager->GetTopVolume()->FindNode("DDIP_1");
  TGeoNode *node2 = gGeoManager->GetTopVolume()->FindNode("YOUT1_1");
  TGeoNode *node3 = gGeoManager->GetTopVolume()->FindNode("YOUT2_1");

  Reve::GeoTopNodeRnrEl* re1 = new Reve::GeoTopNodeRnrEl(gGeoManager,node1);
  re1->UseNodeTrans();
  gReve->AddGlobalRenderElement(re1);

  Reve::GeoTopNodeRnrEl* re2 = new Reve::GeoTopNodeRnrEl(gGeoManager,node2);
  re2->UseNodeTrans();
  gReve->AddGlobalRenderElement(re2);

  Reve::GeoTopNodeRnrEl* re3 = new Reve::GeoTopNodeRnrEl(gGeoManager,node3);
  re3->UseNodeTrans();
  gReve->AddGlobalRenderElement(re3);

  gReve->Redraw3D(kTRUE);

}
