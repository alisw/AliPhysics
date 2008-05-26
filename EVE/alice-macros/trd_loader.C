//
// How to steer the TRD loaders from a macro
// For the usage of only the TRD data containers and 
// AliEve event loop check the macro "trd_detectors.C"
// 
// Usage:
// .L trd_loader.C
// AliEveTRDLoader *loader = trd_loader();
// loader->NextEvent();
// 
// Caution:
// In order to update the screen one has to go to GLViewer 
// and click the UpdateScene button after each NextEvent().
// 
// Author:
// Alex Bercuci (A.Bercuci@gsi.de)
//
AliEveTRDLoader* trd_loader()
{
  // init MC loader
  AliEveTRDLoader *loader = new AliEveTRDLoaderSim("MC");
  
  // init single file loader
  // AliEveTRDLoader *loader = new AliEveTRDLoader("Digits");

  // link the run loader and define the chamber setting and data type
  loader->Open("galice.root");
  loader->AddChambers(3, 3);
  loader->SetDataType(AliEveTRDLoader::kTRDHits | AliEveTRDLoader::kTRDDigits | AliEveTRDLoader::kTRDClusters);

  // load first event
  loader->GoToEvent(0);
  
  // register loader with alieve
  gEve->AddElement(loader);
  loader->SpawnEditor();
  gEve->Redraw3D();

  return loader;
}
