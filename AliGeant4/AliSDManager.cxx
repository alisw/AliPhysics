// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "AliSDManager.h"
//#include "AliSDMessenger.h"
#include "AliSensitiveDetector.h"
#include "AliLegoSensitiveDetector.h"
#include "AliGlobals.h"
#include "AliRun.h"
#include "AliModule.h"

#include "TG4GeometryManager.h"

#include <G4SDManager.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>

AliSDManager* AliSDManager::fgInstance = 0;

AliSDManager* AliSDManager::Instance()
{
// Returns the singleton instance.
// Creates the instance if it does not exist.
// ---

  if (!fgInstance) new AliSDManager();
  
  return fgInstance;
}  
    
AliSDManager::AliSDManager()
  : fNofLVWithSD(0)
{
//
  //fMessenger = new AliSDMessenger(this);
  fgInstance = this;
}

AliSDManager::AliSDManager(const AliSDManager& right) {
//
  AliGlobals::Exception(
    "Singleton AliSDManager is protected from copying.");
}

AliSDManager::~AliSDManager()
{
  //delete fMessenger;
}

// operators

AliSDManager& AliSDManager::operator=(const AliSDManager& right)
{
  // check assignement to self
  if (this == &right) return *this;

  AliGlobals::Exception(
     "Singleton AliSDManager is protected from assigning.");
    
  return *this;  
}    
          
// private methods

void AliSDManager::CreateLegoSD(G4LogicalVolume* lv, AliLego* lego) const
{ 
// Replaces the existing sensitive detector of the logical volume
// with the lego sensitive detector.
// ---

  TG4GeometryManager* pGeometryManager = TG4GeometryManager::Instance();
  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();

  G4String lvName = lv->GetName(); 
  G4String sdName = "/Alice/lego/" + lvName;

  // cut copy number from sdName
  pGeometryManager->G4ToG3VolumeName(sdName);
  
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

void AliSDManager::UnsetLegoSD(G4LogicalVolume* lv) const
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
      G4String text = "AliSDManager::UnsetLego: \n";
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

void AliSDManager::CreateSD(G4LogicalVolume* lv, AliModule* module) const
{ 
// Creates/retrieves a sensitive detector for the logical volume.
// ---

  TG4GeometryManager* pGeometryManager = TG4GeometryManager::Instance();
  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();

  G4String lvName = lv->GetName(); 
  
  G4String moduleName = module->GetName();
  G4String sdName = "/Alice/" + moduleName + "/" + lvName;

  // cut copy number from sdName
  pGeometryManager->G4ToG3VolumeName(sdName);
  
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

AliModule* AliSDManager::FindAliModule(G4LogicalVolume* lv) const
{
// Finds the module containing specified logical volume.
// ---

  // geometry manager
  TG4GeometryManager* pGeometryManager = TG4GeometryManager::Instance();

  // get g3 volume name
  G4String g3Name = lv->GetName();
  pGeometryManager->G4ToG3VolumeName(g3Name);
  
  // get module name from the map
  G4String moduleName = pGeometryManager->GetMapSecond(g3Name);
  
  // find module from gAlice
  AliModule* module = gAlice->GetModule(moduleName);
  if (!module) {
    G4String text = "AliSDManager::FindAliModule:\n";
    text = text + "    AliModule " + moduleName;
    text = text + " (mapped from logical volume " + lv->GetName() + ")\n";
    text = text + "    has not been found in gAlice.";
    AliGlobals::Exception(text);
  }  
   
  return module;
}  

void AliSDManager::SetLego(AliLego* lego) const 
{ 
// Replaces the existing sensitive detectors 
// with lego sensitive detectors.
// ---

  if (!lego) {
    G4String text = "AliSDManager::SetLego: \n";
    text = text + "    No AliLego is defined.";
    AliGlobals::Warning(text);
    return;
  }  

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  G4int nofLV = pLVStore->entries();
  
  for (G4int i=0; i<nofLV; i++) {
    G4LogicalVolume* lv = (*pLVStore)[i];
    CreateLegoSD(lv, lego);
  }
}

void AliSDManager::UnsetLego() const
{
// Replace the lego sensitive detectors 
// back with standard sensitive detectors.
// ---

  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  G4int nofLV = pLVStore->entries();
  
  // set back standard sensitive detectors
  for (G4int i=0; i<nofLV; i++) {
    G4LogicalVolume* lv = (*pLVStore)[i];
    UnsetLegoSD(lv);
  }  

  // The lego sensitive detectors are not deleted here
  // as there is no way how to unregister them
  // from G4SDManager
}
