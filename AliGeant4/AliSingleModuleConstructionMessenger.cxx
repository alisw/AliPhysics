// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "AliSingleModuleConstructionMessenger.h"
#include "AliSingleModuleConstruction.h"
#include "AliGlobals.h"

#include <G4UIcmdWithABool.hh>

AliSingleModuleConstructionMessenger::AliSingleModuleConstructionMessenger(
   AliSingleModuleConstruction* moduleConstruction, G4String moduleName)
 : fModuleConstruction(moduleConstruction)
{
//
  G4String dirName = "/aliDet/"; 
  dirName = dirName + moduleName + "/"; 

  G4String commandPath = dirName + "setAllSensitive";
  fSetAllSensitiveCmd = new G4UIcmdWithABool(commandPath, this);
  G4String guidance =   "If true: set all " + moduleName;
  guidance = guidance + "logical volumes sensitive.\n";
  guidance = guidance + "      (Each logical is volume associated with ";
  guidance = guidance + "a sensitive detector.)";    
  guidance = guidance + "If false: only volumes defined with a sensitive ";
  guidance = guidance + "tracking medium\n";
  guidance = guidance + "are associated with a sensitive detector.";
  fSetAllSensitiveCmd->SetGuidance(guidance);
  fSetAllSensitiveCmd->SetParameterName("detSensitivity", false);
  fSetAllSensitiveCmd->AvailableForStates(PreInit, Init);  
}

AliSingleModuleConstructionMessenger::AliSingleModuleConstructionMessenger() {
//
}

AliSingleModuleConstructionMessenger::AliSingleModuleConstructionMessenger(
                            const AliSingleModuleConstructionMessenger& right) {
//
  AliGlobals::Exception(
    "AliSingleModuleConstructionMessenger is protected from copying.");
}

AliSingleModuleConstructionMessenger::~AliSingleModuleConstructionMessenger() {
//
  delete fSetAllSensitiveCmd;  
}

// operators

AliSingleModuleConstructionMessenger& 
AliSingleModuleConstructionMessenger::operator=(
                            const AliSingleModuleConstructionMessenger& right)
{
  // check assignement to self
  if (this == &right) return *this;

  AliGlobals::Exception(
     "AliSingleModuleConstructionMessenger is protected from assigning.");
    
  return *this;  
}    
          
// public methods
  
void AliSingleModuleConstructionMessenger::SetNewValue(G4UIcommand* command,
                                                 G4String newValues)
{
// Applies command to the associated object.
// ---

  if (command == fSetAllSensitiveCmd) {
    fModuleConstruction
      ->SetAllLVSensitive(fSetAllSensitiveCmd->GetNewBoolValue(newValues)); 
  } 
}

