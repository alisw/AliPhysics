// $Id$

void MUON_geom()
{
  using namespace std;

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  TString fileName = rl->GetFileName();
  Int_t length = fileName.Length();
  fileName.Resize(length-11);
  fileName.Append("geometry.root");

  gGeoManager = gReve->GetGeometry(fileName.Data());

  TGeoNode *node1 = gGeoManager->GetTopVolume()->FindNode("DDIP_1");
  TGeoNode *node2 = gGeoManager->GetTopVolume()->FindNode("YOUT1_1");
  TGeoNode *node3 = gGeoManager->GetTopVolume()->FindNode("YOUT2_1");

  Reve::GeoTopNodeRnrEl* re1 = new Reve::GeoTopNodeRnrEl(gGeoManager,node1);
  re1->SetUseNodeTrans(kTRUE);
  gReve->AddGlobalRenderElement(re1);

  Reve::GeoTopNodeRnrEl* re2 = new Reve::GeoTopNodeRnrEl(gGeoManager,node2);
  re2->SetUseNodeTrans(kTRUE);
  gReve->AddGlobalRenderElement(re2);

  Reve::GeoTopNodeRnrEl* re3 = new Reve::GeoTopNodeRnrEl(gGeoManager,node3);
  re3->SetUseNodeTrans(kTRUE);
  gReve->AddGlobalRenderElement(re3);

  gReve->Redraw3D();

}
