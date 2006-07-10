// $Id$

void geom_phos()
{
  using namespace std;

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");

  for(Int_t i=1; i<=5; ++i) {
    TGeoNode* node = gGeoManager->GetTopVolume()->FindNode(Form("PHOS_%d", i));

    Reve::GeoTopNodeRnrEl* re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
    re->SetUseNodeTrans(kTRUE);
    gReve->AddGlobalRenderElement(re);
  }

  gReve->Redraw3D();
}
