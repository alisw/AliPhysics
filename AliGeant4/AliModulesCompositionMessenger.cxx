// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliModulesCompositionMessenger
// ------------------------------------
// See the class description in the header file.

#include "AliModulesCompositionMessenger.h"
#include "AliModulesComposition.h"
#include "AliGlobals.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>

//_____________________________________________________________________________
AliModulesCompositionMessenger::AliModulesCompositionMessenger(
                                   AliModulesComposition* modulesComposition)
  : fModulesComposition(modulesComposition)
{
//
  fDirectory = new G4UIdirectory("/aliDet/");
  fDirectory->SetGuidance("Detector construction control commands.");

  fFieldTypeCmd = new G4UIcmdWithAString("/aliDet/fieldType", this);
  G4String guidance =   "Select type of magnetic field:\n";
  guidance = guidance + "  MCApplication:  field defined by MC application (default)\n";
  guidance = guidance + "  Uniform:        uniform magnetic field\n";
  guidance = guidance + "  None:           no magnetic field";
  fFieldTypeCmd->SetGuidance(guidance);
  fFieldTypeCmd->SetParameterName("FieldType", true);
  fFieldTypeCmd->SetCandidates("MCApplication Uniform None");   
  fFieldTypeCmd->SetDefaultValue("MCApplication");
  fFieldTypeCmd->AvailableForStates(PreInit);

  fUniformFieldValueCmd 
    = new G4UIcmdWithADoubleAndUnit("/aliDet/uniformFieldValue", this);
  fUniformFieldValueCmd
    ->SetGuidance("Define uniform magnetic field in Z direction.");
  fUniformFieldValueCmd
    ->SetGuidance("(Uniform magnetic field type has to be selected first.)");
  fUniformFieldValueCmd->SetParameterName("UniformFieldValue", false, false);
  fUniformFieldValueCmd->SetDefaultUnit("tesla");
  fUniformFieldValueCmd->SetUnitCategory("Magnetic flux density");
  fUniformFieldValueCmd->AvailableForStates(Idle);  
  
  fSetReadGeometryCmd 
    = new G4UIcmdWithABool("/aliDet/readGeometry", this);
  fSetReadGeometryCmd->SetGuidance("Read geometry from g3calls.dat files");
  fSetReadGeometryCmd->SetParameterName("readGeometry", false);
  fSetReadGeometryCmd->AvailableForStates(PreInit);  
 
  fSetWriteGeometryCmd 
    = new G4UIcmdWithABool("/aliDet/writeGeometry", this);
  fSetWriteGeometryCmd->SetGuidance("Write geometry to g3calls.dat file");
  fSetWriteGeometryCmd->SetParameterName("writeGeometry", false);
  fSetWriteGeometryCmd->AvailableForStates(PreInit);   

  fPrintMaterialsCmd 
    = new G4UIcmdWithoutParameter("/aliDet/printMaterials", this);
  fPrintMaterialsCmd->SetGuidance("Prints all materials.");
  fPrintMaterialsCmd->AvailableForStates(PreInit, Init, Idle);   

  fGenerateXMLCmd 
    = new G4UIcmdWithoutParameter("/aliDet/generateXML", this);
  fGenerateXMLCmd->SetGuidance("Generate geometry XML file.");
  fGenerateXMLCmd->AvailableForStates(Idle);   
}

//_____________________________________________________________________________
AliModulesCompositionMessenger::AliModulesCompositionMessenger() {
//
}

//_____________________________________________________________________________
AliModulesCompositionMessenger::AliModulesCompositionMessenger(
                                const AliModulesCompositionMessenger& right)
{
//
  AliGlobals::Exception(
    "AliModulesCompositionMessenger is protected from copying.");
}

//_____________________________________________________________________________
AliModulesCompositionMessenger::~AliModulesCompositionMessenger() {
//
  delete fDirectory;
  delete fFieldTypeCmd;
  delete fUniformFieldValueCmd;
  delete fSetReadGeometryCmd;
  delete fSetWriteGeometryCmd;
  delete fPrintMaterialsCmd;
  delete fGenerateXMLCmd;
}

// operators

//_____________________________________________________________________________
AliModulesCompositionMessenger& 
AliModulesCompositionMessenger::operator=(
                                const AliModulesCompositionMessenger& right)
{
  // check assignement to self
  if (this == &right) return *this;

  AliGlobals::Exception(
     "AliModulesCompositionMessenger is protected from assigning.");
    
  return *this;  
}    
          
// public methods
  
//_____________________________________________________________________________
void AliModulesCompositionMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
// Applies command to the associated object.
// ---

  if( command == fFieldTypeCmd ) { 
    if (newValues == "MCApplication") 
      fModulesComposition->SetFieldType(kMCApplicationField); 
    if (newValues == "Uniform") 
      fModulesComposition->SetFieldType(kUniformField); 
    if (newValues == "None") 
      fModulesComposition->SetFieldType(kNoField); 
  }
  if (command == fUniformFieldValueCmd) {  
    fModulesComposition
      ->SetUniformFieldValue(fUniformFieldValueCmd->GetNewDoubleValue(newValues)); 
  }
  else if (command == fSetReadGeometryCmd) {
    fModulesComposition->SetReadGeometry(
                         fSetReadGeometryCmd->GetNewBoolValue(newValues));
  }  
  else if (command == fSetWriteGeometryCmd) {
    fModulesComposition->SetWriteGeometry(
                         fSetWriteGeometryCmd->GetNewBoolValue(newValues));
  }    
  else if (command == fPrintMaterialsCmd) {
    fModulesComposition->PrintMaterials();
  }    
  else if (command == fGenerateXMLCmd) {
    fModulesComposition->GenerateXMLGeometry();
  }    
}
