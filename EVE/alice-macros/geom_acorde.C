// $Id$

void geom_acorde()
{
  using namespace std;

  gGeoManager = gReve->GetGeometry("geometry.root");

  Reve::RenderElementList* list = new Reve::RenderElementList("ACORDE");
  gReve->AddGlobalRenderElement(list);

  for(Int_t i=1; i<61; ++i) {
    char form[10000];
    sprintf(form, "ACORDE1_%d", i);
    TGeoNode* node = gGeoManager->GetTopVolume()->FindNode(form);       
    Reve::GeoTopNodeRnrEl* re =  new Reve::GeoTopNodeRnrEl(gGeoManager, node);
    re->UseNodeTrans();
    gReve->AddGlobalRenderElement(list, re);
    // gReve->AddGlobalRenderElement(re, list); // For EVE-dev
  }
  
  gReve->Redraw3D();  
}
