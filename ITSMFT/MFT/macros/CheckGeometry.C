#if !defined(__CINT__) || defined(__MAKECINT__)
#include "ARVersion.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliGeomManager.h"
#include "AliMC.h"
#include <TROOT.h>
#include "AliRun.h"
#include <TGeoManager.h>
#include <TString.h>
#include <TInterpreter.h>
#endif

void CheckGeometry(const char* cfgFile) {

  // Produce the ideal geometry 
  
  AliCDBManager* cdb = AliCDBManager::Instance();
  AliCDBStorage* storage = 0;
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  
  if(!gSystem->AccessPathName("geometry.root")){
    Printf("Deleting existing \"geometry.root\"");
    gSystem->Exec("rm -rf geometry.root");
  }
  gSystem->Load("libgeant321");
  gSystem->Load("libqpythia");
  gSystem->Load("libAliPythia6");
  gROOT->LoadMacro(cfgFile);
  gInterpreter->ProcessLine(gAlice->GetConfigFunction());
  
  gAlice->GetMCApp()->Init();
  
  if(!gGeoManager){
    Printf("Unable to produce a valid geometry to be put in the CDB!");
    return;
  }
  
  if(gSystem->AccessPathName("geometry.root")){
    Printf("Did not find freshly written \"geometry.root\" file. Exiting ...");
    return;
  }
  
  Printf("Reloading freshly written geometry.root file");
  if (TGeoManager::IsLocked()) TGeoManager::UnlockGeometry();
  AliGeomManager::LoadGeometry("geometry.root");

  gGeoManager->DefaultColors(); // assign default colors according to Z of material

  TIter next(gAlice->Modules());
  AliModule *detector;
  AliMFTGeomTGeo *geomTGeo;
  while((detector = dynamic_cast<AliModule*>(next()))) {
    printf(Form("Detector: %s \n", detector->GetName()));
    if (strcmp(detector->GetName(),"MFT") == 0) {
      AliMFT *mft = (AliMFT*)(detector);
      geomTGeo = mft->GetGeomTGeo();
      printf("MFT has: %d disks.\n",geomTGeo->GetNDisks());
    }
  }

}


