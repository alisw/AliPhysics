// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliModulesComposition
// ---------------------------
// See the class description in the header file.

#include "AliModulesComposition.h"
#include "AliModuleConstruction.h"
#include "AliMagneticField.h"
#include "AliGlobals.h"
#include "AliFiles.h"
#include "AliModule.h"

#include "TG4GeometryManager.h"

#include <G4Material.hh>

//_____________________________________________________________________________
AliModulesComposition::AliModulesComposition()
  : AliVerbose("modulesComposition"),
    fReadGeometry(false),
    fWriteGeometry(false),
    fMagneticField(0),
    fMessenger(this) {
//
}

//_____________________________________________________________________________
AliModulesComposition::AliModulesComposition(const AliModulesComposition& right)
  : AliVerbose("modulesComposition"),
    fMessenger(this)
{
//
  AliGlobals::Exception("AliModulesComposition is protected from copying.");  
}

//_____________________________________________________________________________
AliModulesComposition::~AliModulesComposition() {
//   
  delete fMagneticField;
}

// operators

//_____________________________________________________________________________
AliModulesComposition& 
AliModulesComposition::operator=(const AliModulesComposition& right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliModulesComposition is protected from assigning.");  

  return *this;  
}    
          
//
// private methods
//

//_____________________________________________________________________________
void AliModulesComposition::Configure()
{ 
// Executes the detectors setup Root macros
// (extracted from AliRoot Config.C) and
// G4 macros.
// ---
  
  // number of modules
  G4int nofModules = fModuleConstructionVector.size();

  if (nofModules == 0) {
    AliGlobals::Warning(
      "AliModulesComposition::Configure: No modules are defined.");
    return;  
  }
  
  for (G4int i=0; i<nofModules; i++) 
    fModuleConstructionVector[i]->Configure();
}      


//_____________________________________________________________________________
void AliModulesComposition::CreateG4Geometry()
{ 
// Constructs geometry.
// G3 tables are process for all modules alltogether.
// --

  // number of modules
  G4int nofModules = fModuleConstructionVector.size();

  if (nofModules == 0) {
    AliGlobals::Warning(
      "AliModulesComposition::CreateG4Geometry: No modules are defined.");
    return;  
  }

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
	
      G4String text = "AliModulesComposition::Construct - Limitation:\n";
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

/*
  // moved to AliSDConstruction
  // to be performed after Init  (required for MUON)
  
  for (i=0; i<nofModules; i++) {

    // construct geometry for display
    fModuleConstructionVector[i]->GetAliModule()->BuildGeometry();
  }
*/

  // reset TG4GeometryManager 
  pGeometryManager->ClearG3Tables();
}

//_____________________________________________________________________________
void AliModulesComposition::SetReadGeometryToModules(G4bool readGeometry)
{
// Sets readGeometry control to all modules.
// ---

  for (G4int i=0; i<fModuleConstructionVector.size(); i++) 
    fModuleConstructionVector[i]->SetReadGeometry(readGeometry);
}    
  
//_____________________________________________________________________________
void AliModulesComposition::SetWriteGeometryToModules(G4bool writeGeometry)
{
// Sets writeGeometry control to all modules.
// ---

  for (G4int i=0; i<fModuleConstructionVector.size(); i++) 
    fModuleConstructionVector[i]->SetWriteGeometry(writeGeometry);
}    

//
// protected methods
//

//_____________________________________________________________________________
void AliModulesComposition::AddModule(const G4String& name, 
                                      G4int version, 
			              AliModuleType moduleType)
{
// Adds module to the module construction vector.
// ---

  AliModuleConstruction* moduleConstruction 
    = new AliModuleConstruction(name, version, moduleType);

  fModuleConstructionVector.push_back(moduleConstruction);
}  
		  		  
//_____________________________________________________________________________
void AliModulesComposition::ConstructModules()
{
// Construct geometry of all modules (both standalone and dependent.)
// ---

  // set common options
  SetReadGeometryToModules(fReadGeometry);
  SetWriteGeometryToModules(fWriteGeometry);
  
  // configure modules
  Configure();

  // construct dependent modules
  CreateG4Geometry();
}  

//_____________________________________________________________________________
void AliModulesComposition::SetProcessConfigToModules(G4bool processConfig)
{
// Sets processConfig control to all modules.
// ---

  for (G4int i=0; i<fModuleConstructionVector.size(); i++) 
    fModuleConstructionVector[i]->SetProcessConfig(processConfig);
}    

//
// public methods
//

//_____________________________________________________________________________
void AliModulesComposition::PrintMaterials() const
{
// Prints all materials.
// ---

  const G4MaterialTable* matTable = G4Material::GetMaterialTable();
  G4cout << *matTable;
}

//_____________________________________________________________________________
void AliModulesComposition::SetMagField(G4double fieldValue)
{
// Sets uniform magnetic field to specified value.
// ---

  // create fields if it does not exist
  if (!fMagneticField) fMagneticField = new AliMagneticField();
  
  // set value
  fMagneticField->SetFieldValue(fieldValue);
}

