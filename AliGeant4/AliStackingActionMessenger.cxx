// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliStackingActionMessenger.h"
#include "AliStackingAction.h"
#include "AliGlobals.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithoutParameter.hh>

AliStackingActionMessenger::AliStackingActionMessenger(
                              AliStackingAction* stackingAction)
  :fStackingAction(stackingAction)
{
// 
  fStackingDirectory = new G4UIdirectory("/aliStacking/");
  fStackingDirectory->SetGuidance("AliStackingAction control commands.");

  fClearStackCmd = new G4UIcmdWithoutParameter("/aliStacking/clearStack", this);
  fClearStackCmd->SetGuidance("Clears the primary stack.");
  fClearStackCmd->AvailableForStates(EventProc);

  fVerboseCmd = new G4UIcmdWithAnInteger("/aliStacking/verbose", this);
  fVerboseCmd->SetGuidance("Set verbose level for AliStackingAction");
  fVerboseCmd->SetParameterName("VerboseLevel", true);
  fVerboseCmd->SetDefaultValue(0);
  fVerboseCmd->SetRange("VerboseLevel >= 0 && VerboseLevel <= 1");
  fVerboseCmd->AvailableForStates(Idle);
}

AliStackingActionMessenger::AliStackingActionMessenger() {
//
}

AliStackingActionMessenger::AliStackingActionMessenger(
                                 const AliStackingActionMessenger& right) {
//				 
  AliGlobals::Exception("AliStackingActionMessenger is protected from copying.");
}

AliStackingActionMessenger::~AliStackingActionMessenger() {
//
  delete fStackingDirectory;
  delete fClearStackCmd;
  delete fVerboseCmd;
}

// operators

AliStackingActionMessenger& 
AliStackingActionMessenger::operator=(const AliStackingActionMessenger &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception(
    "AliStackingActionMessenger is protected from assigning.");

  return *this;
}

// public methods

void AliStackingActionMessenger::SetNewValue(G4UIcommand* command, 
       G4String newValue)
{ 
// Applies command to the associated object.
// ---

  if (command == fClearStackCmd) { 
    fStackingAction->ClearPrimaryStack(); 
  }
  else if (command == fVerboseCmd) { 
    fStackingAction
      ->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue)); 
  }  
}
