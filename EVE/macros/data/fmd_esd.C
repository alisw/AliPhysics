
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TEveManager.h>

#include <AliEveFMDLoader.h>
#endif

void fmd_esd()
{
  gStyle->SetPalette(1);
  AliEveFMDLoader* gFmdLoader = AliEveFMDLoader::Instance();
  gFmdLoader->LoadESD();
  gEve->Redraw3D();
}
