void fmd_digits()
{
  gStyle->SetPalette(1);
  AliEveFMDLoader* gFmdLoader = AliEveFMDLoader::Instance();
  gFmdLoader->LoadDigits();
  gEve->Redraw3D();
}
