
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TEveManager.h>

#include <AliEveFMDLoader.h>
#endif

void fmd_raw()
{
    printf("*** RAW FMD ***");
    
  gStyle->SetPalette(1);
  AliEveFMDLoader* gFmdLoader = AliEveFMDLoader::Instance();
  gFmdLoader->LoadRaw();
  gEve->Redraw3D();
}
