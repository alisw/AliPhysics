// $Id$

void MUON_geomAll()
{

  using namespace std;

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  TString fileName = rl->GetFileName();
  Int_t length = fileName.Length();
  fileName.Resize(length-11);
  fileName.Append("geometry.root");

  gGeoManager = gReve->GetGeometry(fileName.Data());

  Reve::GeoTopNodeRnrEl* topn_re = new Reve::GeoTopNodeRnrEl
    (gGeoManager, gGeoManager->GetTopNode());
  gReve->AddGlobalRenderElement(topn_re);
  gReve->DrawRenderElement(topn_re);

}
