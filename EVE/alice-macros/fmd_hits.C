void fmd_hits()
{
  gStyle->SetPalette(1);
  AliEveFMDLoader* gFmdLoader = AliEveFMDLoader::Instance();
  gFmdLoader->LoadHits();
  gEve->Redraw3D();
}
