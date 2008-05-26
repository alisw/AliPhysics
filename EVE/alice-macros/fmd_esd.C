void fmd_esd()
{
  gStyle->SetPalette(1);
  AliEveFMDLoader* gFmdLoader = AliEveFMDLoader::Instance();
  gFmdLoader->LoadESD();
  gEve->Redraw3D();
}
