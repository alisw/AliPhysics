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

//_____________________________________________________________________________
AliStackingActionMessenger::AliStackingActionMessenger(
                              AliStackingAction* stackingAction)
  :fStackingAction(stackingAction)
{
// 
  fStackingDirectory = new G4UIdirectory("/aliStacking/");
  fStackingDirectory->SetGuidance("AliStackingAction control commands.");

  fVerboseCmd = new G4UIcmdWithAnInteger("/aliStacking/verbose", this);
  fVerboseCmd->SetGuidance("Set verbose level for AliStackingAction");
  fVerboseCmd->SetParameterName("VerboseLevel", true);
  fVerboseCmd->SetDefaultValue(0);
  fVerboseCmd->SetRange("VerboseLevel >= 0 && VerboseLevel <= 1");
  fVerboseCmd->AvailableForStates(Idle);
}

//_____________________________________________________________________________
AliStackingActionMessenger::AliStackingActionMessenger() {
//
}

//_____________________________________________________________________________
AliStackingActionMessenger::AliStackingActionMessenger(
                                 const AliStackingActionMessenger& right) {
//				 
  AliGlobals::Exception("AliStackingActionMessenger is protected from copying.");
}

//_____________________________________________________________________________
AliStackingActionMessenger::~AliStackingActionMessenger() {
//
  delete fStackingDirectory;
  delete fVerboseCmd;
}

// operators

//_____________________________________________________________________________
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

//_____________________________________________________________________________
void AliStackingActionMessenger::SetNewValue(G4UIcommand* command, 
       G4String newValue)
{ 
// Applies command to the associated object.
// ---

  if (command == fVerboseCmd) { 
    fStackingAction
      ->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue)); 
  }  
}
