// $Id$

void geom_emcal()
{
  gGeoManager = gEve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("XEN1_1");

  TEveGeoTopNode* emcal_re = new TEveGeoTopNode(gGeoManager, node);
  gEve->AddGlobalElement(emcal_re);
  gEve->Redraw3D();
}
