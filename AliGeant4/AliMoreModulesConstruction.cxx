// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliMoreModulesConstruction
// --------------------------------
// See the class description in the header file.

#include "AliMoreModulesConstruction.h"
#include "AliSingleModuleConstruction.h"
#include "AliModule.h"
#include "AliGlobals.h"
#include "AliFiles.h"

#include "TG4GeometryManager.h"

#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>

//_____________________________________________________________________________
AliMoreModulesConstruction::AliMoreModulesConstruction() {
//
}

//_____________________________________________________________________________
AliMoreModulesConstruction::AliMoreModulesConstruction(
                               const AliMoreModulesConstruction& right)
{
  // copy stuff
  *this = right;
}  
  			       

//_____________________________________________________________________________
AliMoreModulesConstruction::~AliMoreModulesConstruction()
{
  // delete module constructions
  fModuleConstructionVector.erase(
    fModuleConstructionVector.begin(), fModuleConstructionVector.end());
}

// operators

//_____________________________________________________________________________
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

  return *this;  
}    
          
// public methods

//_____________________________________________________________________________
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
  					   
//_____________________________________________________________________________
void AliMoreModulesConstruction::Configure(const AliFiles& files)
{ 
// Executes the detectors setup Root macros
// (extracted from AliRoot Config.C) and
// G4 macros.
// ---
  
  // number of modules
  G4int nofModules = fModuleConstructionVector.size();

  if (nofModules == 0) {
    AliGlobals::Warning(
      "AliMoreModulesConstruction::Construct(): No modules are defined.");
  }
  else 
    for (G4int i=0; i<nofModules; i++) 
      fModuleConstructionVector[i]->Configure(files);
}      

//_____________________________________________________________________________
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

    G4int i;
    for (i=0; i<nofModules; i++) {

      // fModuleConstructionVector[i]->Configure(files);
    
      // register module name in the name map
      AliModule* module = fModuleConstructionVector[i]->GetAliModule();
      pGeometryManager->SetMapSecond(module->GetName());	

      G4bool readGeometry = fModuleConstructionVector[i]->GetReadGeometry();
      G4bool writeGeometry = fModuleConstructionVector[i]->GetWriteGeometry();
      G4String dataFilePath = fModuleConstructionVector[i]->GetDataFilePath();

      if (readGeometry) {
        // TG4GeometryManager uses g3tog4 methods for reading
	// g3calls.dat files - as these methods do not fill name map
	// they cannot be used for constructing more modules
	// together
	//
        // pGeometryManager->SetWriteGeometry(false);
        // pGeometryManager->ReadG3Geometry(dataFilePath);
	
	G4String text = "AliMoreModulesConstruction::Construct - Limitation:\n";
	text = text + "    Reading g3calls.dat is not implemented.";
	AliGlobals::Exception(text);
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
        
        if (writeGeometry) 
          pGeometryManager->CloseOutFile();
      }	
    }  
  
    // construct G4 geometry
    pGeometryManager->CreateG4Geometry();

    // print name map
    // pGeometryManager->PrintNameMap();

    for (i=0; i<nofModules; i++) {
      // set the detector frame (envelope)
      // (without warning output if enevelope is not defined)
      fModuleConstructionVector[i]->SetDetFrame(false);

      // construct geometry for display
      fModuleConstructionVector[i]->GetAliModule()->BuildGeometry();
    }

    // reset TG4GeometryManager 
    pGeometryManager->ClearG3Tables();

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
