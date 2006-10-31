// $Id$

void geom_rich()
{
  using namespace std;

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");

  Reve::RenderElementList* list = new Reve::RenderElementList("RICH");
  gReve->AddGlobalRenderElement(list);

  for(Int_t i=1; i<=7; ++i) {
    TGeoNode* node = gGeoManager->GetTopVolume()->FindNode(Form("RICH_%d", i));

    Reve::GeoTopNodeRnrEl* re = new Reve::GeoTopNodeRnrEl(gGeoManager, node);
    re->SetUseNodeTrans(kTRUE);
    gReve->AddGlobalRenderElement(list, re);
  }

  gReve->Redraw3D();
}
