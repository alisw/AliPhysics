void DrawGeom( Bool_t checkOverlaps = kFALSE)
{
  TGeoManager::Import("geometry.root");
//  gGeoManager->DefaultColors();

  gGeoManager->GetVolume("ALIC")->InvisibleAll();
  gGeoManager->SetVisLevel(8);
  new TBrowser;
  if(CheckOverlaps){
    gGeoManager->CheckOverlaps(0.1);
    gGeoManager->PrintOverlaps();
  }
  gGeoManager->GetVolume("MFT")->Draw("ogl");
}
