
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TEveManager.h>

#include <AliEveFMDLoader.h>
#endif

void fmd_hits()
{
  gStyle->SetPalette(1);
  AliEveFMDLoader* gFmdLoader = AliEveFMDLoader::Instance();
  gFmdLoader->LoadHits();
  gEve->Redraw3D();
}
