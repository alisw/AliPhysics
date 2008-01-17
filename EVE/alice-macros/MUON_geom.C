// $Id$

void MUON_geom()
{
  TString dataPath = TString(Alieve::gEvent->GetTitle());
  dataPath.Append("/geometry.root");
  gGeoManager = gEve->GetGeometry(dataPath.Data());

  TGeoNode *node1 = gGeoManager->GetTopVolume()->FindNode("DDIP_1");
  TGeoNode *node2 = gGeoManager->GetTopVolume()->FindNode("YOUT1_1");
  TGeoNode *node3 = gGeoManager->GetTopVolume()->FindNode("YOUT2_1");

  TEveGeoTopNode* re1 = new TEveGeoTopNode(gGeoManager,node1);
  re1->UseNodeTrans();
  gEve->AddGlobalElement(re1);

  TEveGeoTopNode* re2 = new TEveGeoTopNode(gGeoManager,node2);
  re2->UseNodeTrans();
  gEve->AddGlobalElement(re2);

  TEveGeoTopNode* re3 = new TEveGeoTopNode(gGeoManager,node3);
  re3->UseNodeTrans();
  gEve->AddGlobalElement(re3);

  gEve->Redraw3D(kTRUE);
}
