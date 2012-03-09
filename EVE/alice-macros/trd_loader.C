//
// How to steer the TRD loaders from a macro
// For the usage of only the TRD data containers and 
// AliEve event loop check the macro "trd_detectors.C"
// 
// Usage:
// .L trd_loader.C
// AliCDBManager *cdb = AliCDBManager::Instance();
// cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
// cdb->SetRun(0)
// AliEveTRDLoader *loader = trd_loader();
// loader->NextEvent();
// loop(loader)
// 
// Author:
// Alex Bercuci (A.Bercuci@gsi.de)
//

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TString.h>
#include <TSystem.h>
#include <TEveManager.h>
#include <TEveTreeTools.h>

#include <AliEveTRDLoader.h>
#endif

AliEveTRDLoader* trd_loader(Int_t event=0)
{
  // init single file loader
  AliEveTRDLoader *loader = new AliEveTRDLoader("Clusters");

  // link the run loader and define the chamber setting and data type
  loader->Open("TRD.RecPoints.root");
  loader->AddChambers(0);
  loader->AddChambers(8);
  loader->AddChambers(9);
  loader->AddChambers(17);
  loader->SetDataType(AliEveTRDLoader::kTRDClusters);

  // load first event
  loader->GoToEvent(event);
  
  // register loader with alieve
  gEve->AddElement(loader);
  loader->SpawnEditor();
  gEve->Redraw3D();

  return loader;
}


void loop(AliEveTRDLoader *loader)
{
  while(loader->NextEvent()){ 
    printf("Event[%d]\n", loader->GetEvent());
    gEve->Redraw3D();
    gSystem->ProcessEvents();
    gSystem->Sleep(2000);
  }
}
