// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "AliMoreModulesConstruction.h"
#include "AliSingleModuleConstruction.h"
#include "AliSDManager.h"
#include "AliSensitiveDetector.h"
#include "AliModule.h"
#include "AliRun.h"
#include "AliGlobals.h"

#include "TG4GeometryManager.h"

#include <G4SDManager.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>

#include <TROOT.h> 
#include <TCint.h> 

AliMoreModulesConstruction::AliMoreModulesConstruction() {
//
  fSDManager = AliSDManager::Instance();
}

AliMoreModulesConstruction::AliMoreModulesConstruction(
                               const AliMoreModulesConstruction& right)
{
  // delete current module constructions
  fModuleConstructionVector.erase(
    fModuleConstructionVector.begin(), fModuleConstructionVector.end());
    
  // create new module constructions   
  G4int nofModules = right.fModuleConstructionVector.size();
  for (G4int i=0; i<nofModules; i++) {
    G4String name = right.fModuleConstructionVector[i]->GetDetName();
    G4int version = right.fModuleConstructionVector[i]->GetVersion();
    AliModuleType type = right.fModuleConstructionVector[i]->GetType();
    AddModule(name, version, type);
  }  
  
  fSDManager = right.fSDManager;
}  
  			       

AliMoreModulesConstruction::~AliMoreModulesConstruction()
{
  // delete module constructions
  fModuleConstructionVector.erase(
    fModuleConstructionVector.begin(), fModuleConstructionVector.end());
}

// operators

AliMoreModulesConstruction& 
AliMoreModulesConstruction::operator=(const AliMoreModulesConstruction& right)
{
  // check assignement to self
  if (this == &right) return *this;

  // delete current module constructions
  fModuleConstructionVector.erase(
    fModuleConstructionVector.begin(), fModuleConstructionVector.end());
    
  // create new module constructions   
  G4int nofModules = right.fModuleConstructionVector.size();
  for (G4int i=0; i<nofModules; i++) {
    G4String name = right.fModuleConstructionVector[i]->GetDetName();
    G4int version = right.fModuleConstructionVector[i]->GetVersion();
    AliModuleType type = right.fModuleConstructionVector[i]->GetType();
    AddModule(name, version, type);
  }  
  
  fSDManager = right.fSDManager;

  return *this;  
}    
          
// private methods

void AliMoreModulesConstruction::CreateSensitiveDetectors(
                                   G4bool allLVSensitive)
{
// Creates sensitive detectors.
// ---

  if (allLVSensitive)
    CreateSensitiveDetectors1();
  else
    CreateSensitiveDetectors2();

  // set static number of logical volumes already processed
  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  fSDManager->SetNofLVWithSD(pLVStore->entries());  
}    

void AliMoreModulesConstruction::CreateSensitiveDetectors1()
{ 
// Creates sensitive detectors.
// Sensitive detectors are set to all logical volumes
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  G4int nofLV = pLVStore->entries();
  
  G4int nofLVWithSD = fSDManager->GetNofLVWithSD();
  for (G4int i=nofLVWithSD; i<nofLV; i++) {
    G4LogicalVolume* lv = (*pLVStore)[i];
    AliModule* module = fSDManager->FindAliModule(lv);
    fSDManager->CreateSD(lv, module);
  }
}

void AliMoreModulesConstruction::CreateSensitiveDetectors2()
{ 
// Creates sensitive detectors.
// Sensitive detectors are set only to logical volumes
// in G3SensVolVector.
// ---

  TG4GeometryManager* pGeometryManager = TG4GeometryManager::Instance();

  G3SensVolVector pSVVector
    = pGeometryManager->GetG3SensVolVector();

  G4int nofSV = pSVVector.entries();
  if (nofSV>0)
    for (G4int isv=0; isv<nofSV; isv++) {
      G4LogicalVolume* lv = pSVVector[isv];
      AliModule* module = fSDManager->FindAliModule(lv);    
      fSDManager->CreateSD(lv, module);
    } 
}

// public methods

void AliMoreModulesConstruction::AddModule(G4String moduleName, G4int version,
                                           AliModuleType moduleType)
{					   
// Adds module specified by name, version and type.
// ---

  // create module construction
  AliSingleModuleConstruction* moduleConstruction 
    = new AliSingleModuleConstruction(moduleName, version, moduleType);

  // add module, module construction to vectors
  fModuleConstructionVector.push_back(moduleConstruction);
}  
  					   
void AliMoreModulesConstruction::Construct()
{ 
// Constructs geometry.
// G3 tables are process for all modules alltogether.
// --

  // number of modules
  G4int nofModules = fModuleConstructionVector.size();

  if (nofModules == 0) {
    AliGlobals::Warning(
      "AliMoreModulesConstruction::Construct(): No modules are defined.");
  }
  else {      
    // get geometry manager
    TG4GeometryManager* pGeometryManager = TG4GeometryManager::Instance();

    G4bool allLVSensitive = false;    
    G4int i;
    for (i=0; i<nofModules; i++) {

      fModuleConstructionVector[i]->Configure();
    
      // register module name in the name map
      AliModule* module = fModuleConstructionVector[i]->GetAliModule();
      pGeometryManager->SetMapSecond(module->GetName());	

      G4bool readGeometry = fModuleConstructionVector[i]->GetReadGeometry();
      G4bool writeGeometry = fModuleConstructionVector[i]->GetWriteGeometry();
      G4String dataFilePath = fModuleConstructionVector[i]->GetDataFilePath();

      if (readGeometry) {
        // create G3 geometry from g3calls.dat
        pGeometryManager->SetWriteGeometry(false);
        pGeometryManager->ReadG3Geometry(dataFilePath);
      }
      else {	            
        // set geometry output stream for this module
        pGeometryManager->SetWriteGeometry(writeGeometry);
        if (writeGeometry) 
          pGeometryManager->OpenOutFile(dataFilePath);
        
        // create geometry from AliRoot
      
        // construct materials
        module->CreateMaterials();   

        // construct G3 geometry
        module->CreateGeometry();
      }	

      // all logical volumes will be made sensitive if any
      // module requires this
      if (fModuleConstructionVector[i]->GetAllLVSensitive()) 
         allLVSensitive = true;
    }  
  
    // construct G4 geometry
    pGeometryManager->CreateG4Geometry();

    // print name map
    // pGeometryManager->PrintNameMap();
    
    // create sensitive detectors
    CreateSensitiveDetectors(allLVSensitive);
  
    for (i=0; i<nofModules; i++) {
      // set the detector frame (envelope)
      // (without warning output if enevelope is not defined)
      fModuleConstructionVector[i]->SetDetFrame(false);

      // build sensitive detectors table
      fModuleConstructionVector[i]->GetAliModule()->Init();
    }  

    // reset TG4GeometryManager 
    pGeometryManager->ClearG3Tables();
  
    // print current total number of logical volumes
    G4cout << "Current total number of sensitive volumes: "
           << pGeometryManager->NofVolumes() << endl;

#ifdef ALICE_VISUALIZE
    // set visualization attributes
    for (i=0; i<nofModules; i++) {
      if (fModuleConstructionVector[i]->GetDetFrame()) {
        fModuleConstructionVector[i]->SetDetVisibility(true);
        fModuleConstructionVector[i]->SetDetColour("Yellow");
      }
    }  
#endif
  }
}
