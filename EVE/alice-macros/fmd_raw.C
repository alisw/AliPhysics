void fmd_raw()
{
  gStyle->SetPalette(1);
  AliEveFMDLoader* gFmdLoader = AliEveFMDLoader::Instance();
  gFmdLoader->LoadRaw();
  gEve->Redraw3D();
}
