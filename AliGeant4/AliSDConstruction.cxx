// $Id$
// Category: digits+hits
//
// See the class description in the header file.

#include "AliSDConstruction.h"
#include "AliSensitiveDetector.h"
#include "AliLegoSensitiveDetector.h"
#include "AliGlobals.h"
#include "AliRun.h"
#include "AliModule.h"

#include "TG4GeometryServices.h"

#include <G4SDManager.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>

#include <TObjArray.h>

AliSDConstruction::AliSDConstruction()
  : TG4VSDConstruction() {
//
}

AliSDConstruction::~AliSDConstruction() {
//
}

// private methods

void AliSDConstruction::InitializeModules()
{
// Initializes AliModules.
// ---
  
  TObjArray* modules = gAlice->Modules();
  TIter next(modules);
  AliModule* module;
  while((module = (AliModule*)next())) 
    module->Init();
}  

AliModule* AliSDConstruction::FindAliModule(G4LogicalVolume* lv) const
{
// Finds the module containing specified logical volume.
// ---

  // geometry manager
  TG4GeometryServices* geometryServices = TG4GeometryServices::Instance();

  // get g3 volume name
  G4String g3Name = lv->GetName();
  geometryServices->G4ToG3VolumeName(g3Name);
  
  // get module name from the map
  G4String moduleName = geometryServices->GetMapSecond(g3Name);
  
  // find module from gAlice
  AliModule* module = gAlice->GetModule(moduleName);
  if (!module) {
    G4String text = "AliSDConstruction::FindAliModule:\n";
    text = text + "    AliModule " + moduleName;
    text = text + " (mapped from logical volume " + lv->GetName() + ")\n";
    text = text + "    has not been found in gAlice.";
    AliGlobals::Exception(text);
  }  
   
  return module;
}  

void AliSDConstruction::CreateSD(G4LogicalVolume* lv, AliModule* module) const
{ 
// Creates/retrieves a sensitive detector for the logical volume.
// ---

  TG4GeometryServices* geometryServices = TG4GeometryServices::Instance();
  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();

  G4String lvName = lv->GetName(); 
  
  G4String moduleName = module->GetName();
  G4String sdName = "/Alice/" + moduleName + "/" + lvName;

  // cut copy number from sdName
  geometryServices->G4ToG3VolumeName(sdName);
  
  // create/retrieve the sensitive detector
  G4VSensitiveDetector* sd = 0; 
  sd = pSDManager->FindSensitiveDetector(sdName);
  if (!sd) {
    AliSensitiveDetector* asd 
      = new AliSensitiveDetector(sdName, module);	
    pSDManager->AddNewDetector(asd);
    // add verbose
    G4cout << "Sensitive detector " << sdName << "(" 
           << asd->GetID() << ") has been created." << G4endl;
    sd = asd;  
  }	
  lv->SetSensitiveDetector(sd);	     
}

void AliSDConstruction::CreateLegoSD(G4LogicalVolume* lv, AliLego* lego) const
{ 
// Replaces the existing sensitive detector of the logical volume
// with the lego sensitive detector.
// ---

  TG4GeometryServices* geometryServices = TG4GeometryServices::Instance();
  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();

  G4String lvName = lv->GetName(); 
  G4String sdName = "/Alice/lego/" + lvName;

  // cut copy number from sdName
  geometryServices->G4ToG3VolumeName(sdName);
  
  // retrieve the standard sensitive detector
  G4VSensitiveDetector* sd = lv->GetSensitiveDetector();

  // create/retrieve the lego sensitive detector
  G4VSensitiveDetector* legoVSD = 0; 
  legoVSD = pSDManager->FindSensitiveDetector(sdName);

  if (!legoVSD) {
    AliLegoSensitiveDetector* legoSD
      = new AliLegoSensitiveDetector(sdName, lego, sd);	
    pSDManager->AddNewDetector(legoSD);
    // add verbose
    G4cout << "Lego sensitive detector " << sdName 
         << ") has been created." << G4endl;
    legoVSD = legoSD;  
  }	  
  lv->SetSensitiveDetector(legoVSD);	     
}

void AliSDConstruction::UnsetLegoSD(G4LogicalVolume* lv) const
{ 
// Replace the lego sensitive detector of the logical volume
// with the standard sensitive detector.
// ---

  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();

  // get the lego sensitive detector
  G4VSensitiveDetector* sd = lv->GetSensitiveDetector();
  AliLegoSensitiveDetector* legoSD = 0; 
  if (sd) {
    // check type
    legoSD = dynamic_cast<AliLegoSensitiveDetector*>(sd);
    if (!legoSD) {
      G4String text = "AliSDConstruction::UnsetLego: \n";
      text = text + "    Wrong sensitive detector type.";
      AliGlobals::Exception(text);
    }	
  }
 
  // get the standard sensitive detector
  G4VSensitiveDetector* standardSD = legoSD->GetStandardSD();

  // set the standard sensitive detector
  lv->SetSensitiveDetector(standardSD);	     
}


// public methods

void AliSDConstruction::Construct()
{ 
// Creates sensitive detectors and initialize AliModules.
// Sensitive detectors are set to all logical volumes
// ---
  
  G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
  G4int nofLV = lvStore->size();
  
  for (G4int i=0; i<nofLV; i++) {
    G4LogicalVolume* lv = (*lvStore)[i];
    AliModule* module = FindAliModule(lv);
    CreateSD(lv, module);
  }
  
  InitializeModules();
}

/*
//#include "TG4GeometryManager.h"
//#include <G3SensVolVector.hh>

void AliSDConstruction::CreateSensitiveDetectors2()
{ 
// Creates sensitive detectors.
// Sensitive detectors are set only to logical volumes
// in G3SensVolVector.
// ---

  G3SensVolVector svVector
    = TG4GeometryManager::Instance()->GetG3SensVolVector();

  G4int nofSV = svVector.entries();
  if (nofSV>0)
    for (G4int isv=0; isv<nofSV; isv++) {
      G4LogicalVolume* lv = svVector[isv];
      AliModule* module = FindAliModule(lv);    
      CreateSD(lv, module);
    } 
}
*/

void AliSDConstruction::SetLego(AliLego* lego) const 
{ 
// Replaces the existing sensitive detectors 
// with lego sensitive detectors.
// ---

  if (!lego) {
    G4String text = "AliSDConstruction::SetLego: \n";
    text = text + "    No AliLego is defined.";
    AliGlobals::Warning(text);
    return;
  }  

  G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
  G4int nofLV = lvStore->size();
  
  for (G4int i=0; i<nofLV; i++) {
    G4LogicalVolume* lv = (*lvStore)[i];
    CreateLegoSD(lv, lego);
  }
}

void AliSDConstruction::UnsetLego() const
{
// Replace the lego sensitive detectors 
// back with standard sensitive detectors.
// ---

  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
  G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
  G4int nofLV = lvStore->size();
  
  // set back standard sensitive detectors
  for (G4int i=0; i<nofLV; i++) {
    G4LogicalVolume* lv = (*lvStore)[i];
    UnsetLegoSD(lv);
  }  

  // The lego sensitive detectors are not deleted here
  // as there is no way how to unregister them
  // from G4SDManager
}
