// $Id$

void MUON_geomAll()
{

  using namespace std;

  TString dataPath = TString(Alieve::gEvent->GetTitle());
  dataPath.Append("/geometry.root");
  gGeoManager = gEve->GetGeometry(dataPath.Data());

  TEveGeoTopNode* topn_re = new TEveGeoTopNode
    (gGeoManager, gGeoManager->GetTopNode());
  
  gEve->AddGlobalElement(topn_re);

  gEve->Redraw3D(kTRUE);
  
}
