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
#include "AliGlobals.h"
#include "AliFiles.h"
#include "AliModule.h"

#include "TG4GeometryManager.h"
#include "TG4MagneticField.h"
#include "TG4UniformMagneticField.h"

#include <G4Material.hh>

//_____________________________________________________________________________
AliModulesComposition::AliModulesComposition()
  : AliVerbose("modulesComposition"),
    fMessenger(this),
    fMagneticFieldType(kMCApplicationField), 
    fMagneticField(0),
    fReadGeometry(false),
    fWriteGeometry(false) {
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
void AliModulesComposition::CreateMagneticField()
{
// Creates standard magnetic field (defined by TVirtualMCApplication).
// ---

  switch (fMagneticFieldType) {
  
    case kMCApplicationField:
      fMagneticField = new TG4MagneticField();
      G4cout << "kMCApplicationField" << G4endl;
      break;

    case kUniformField:
      fMagneticField = new TG4UniformMagneticField();
      G4cout << "kUniformField" << G4endl;
      break;
      
    case kNoField:
      G4cout << "kNoField" << G4endl;
      ;;
  }  
}

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

  for (G4int i=0; i<G4int(fModuleConstructionVector.size()); i++) 
    fModuleConstructionVector[i]->SetReadGeometry(readGeometry);
}    
  
//_____________________________________________________________________________
void AliModulesComposition::SetWriteGeometryToModules(G4bool writeGeometry)
{
// Sets writeGeometry control to all modules.
// ---

  for (G4int i=0; i<G4int(fModuleConstructionVector.size()); i++) 
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

  // create magnetic field
  CreateMagneticField();  

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

  for (G4int i=0; i<G4int(fModuleConstructionVector.size()); i++) 
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
void AliModulesComposition::SetFieldType(TG4MagneticFieldType fieldType)
{
// Selects from available magnetic field types:
// field defined by TVirtualMCApplication, uniform magnetic field
// or no magnetic field.
// Applicable only when no field has been yet created (PreInit).
// ---

  if (fMagneticField) {
     G4String text = "AliModulesComposition::SetFieldType :\n";
     text = text + "    The magnetic field already exists.";
     text = text + "    Selection was ignored.";
     AliGlobals::Warning(text);
  }   

  fMagneticFieldType = fieldType;
}  

//_____________________________________________________________________________
void AliModulesComposition::SetUniformFieldValue(G4double fieldValue)
{
// Sets uniform magnetic field to specified value.
// ---
  
  if (!fMagneticField) {
     G4String text = "AliModulesComposition::SetUniformMagField: \n";
     text = text + "    Magnetic field is not defined.";
     AliGlobals::Exception(text);
  }   

  // Check field type 
  TG4UniformMagneticField* uniformField 
    =dynamic_cast<TG4UniformMagneticField*>(fMagneticField);

  if (!uniformField) {
     G4String text = "AliModulesComposition::SetUniformMagField: \n";
     text = text + "    Defined magnetic field is not uniform.";
     AliGlobals::Exception(text);
  }   

  // Set value
  uniformField->SetFieldValue(fieldValue);
}
