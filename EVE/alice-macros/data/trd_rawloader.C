//
// How to steer the TRD loaders from a macro
// For the usage of only the TRD data containers and 
// AliEve event loop check the macro "trd_detectors.C"
// 
// Usage:
// .L trd_rawloader.C
// AliEveTRDLoader *raw = trd_rawloader(filename);
// raw->NextEvent();
// 
// Caution:
// In order to update the screen one has to go to GLViewer 
// and click the UpdateScene button after each NextEvent().
// 
// Author:
// Alex Bercuci (A.Bercuci@gsi.de)
// Minjug Kweon (minjung@physi.uni-heidelberg.de)
//

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TEveManager.h>

#include <AliEveTRDLoader.h>
#include <AliEveTRDLoaderImp.h>
#endif


AliEveTRDLoader* trd_rawloader(Char_t *file)
{
  Int_t fSuperModule = 0; // -1 for all
  Int_t fStack = 4;       // -1 for all
  Int_t fLayer = -1;      // -1 for all

  // init RAW loader
  AliEveTRDLoaderRaw *raw = new AliEveTRDLoaderRaw("RAW");
  raw->SetDataType(AliEveTRDLoader::kTRDRawRoot);
  raw->AddChambers(fSuperModule, fStack, fLayer);
  raw->Open(file);

  // load first event
  raw->GoToEvent(0);
  
  // register raw with alieve
  gEve->AddElement(raw);
  raw->SpawnEditor();
  gEve->Redraw3D();

  return raw;
}
