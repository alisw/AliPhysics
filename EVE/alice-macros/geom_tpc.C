// $Id$

void geom_tpc()
{
  gGeoManager = gEve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("TPC_M_1");

  TEveGeoTopNode* tpc_re = new TEveGeoTopNode(gGeoManager, node);
  gEve->AddGlobalElement(tpc_re);
  gEve->Redraw3D();
}
