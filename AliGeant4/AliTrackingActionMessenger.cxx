// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliTrackingActionMessenger.h"
#include "AliTrackingAction.h"
#include "AliGlobals.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAnInteger.hh>

AliTrackingActionMessenger::AliTrackingActionMessenger(
                              AliTrackingAction* trackingAction)
  :fTrackingAction(trackingAction)
{
// 
  fTrackingDirectory = new G4UIdirectory("/aliTracking/");
  fTrackingDirectory->SetGuidance("AliTrackingAction control commands.");

  fVerboseCmd = new G4UIcmdWithAnInteger("/aliTracking/verbose", this);
  fVerboseCmd->SetGuidance("Set verbose level for AliTrackingAction");
  fVerboseCmd->SetParameterName("VerboseLevel", true);
  fVerboseCmd->SetDefaultValue(2);
  fVerboseCmd->SetRange("VerboseLevel >= 0 && VerboseLevel <= 3");
  fVerboseCmd->AvailableForStates(Idle);
}

AliTrackingActionMessenger::AliTrackingActionMessenger() {
//
}

AliTrackingActionMessenger::AliTrackingActionMessenger(
                                 const AliTrackingActionMessenger& right) {
//				 
  AliGlobals::Exception(
    "AliTrackingActionMessenger is protected from copying.");
}

AliTrackingActionMessenger::~AliTrackingActionMessenger() {
//
  delete fTrackingDirectory;
  delete fVerboseCmd;
}

// operators

AliTrackingActionMessenger& 
AliTrackingActionMessenger::operator=(const AliTrackingActionMessenger &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception(
    "AliTrackingActionMessenger is protected from assigning.");

  return *this;
}

// public methods

void AliTrackingActionMessenger::SetNewValue(G4UIcommand* command, 
       G4String newValue)
{ 
// Applies command to the associated object.
// ---

  if(command == fVerboseCmd)
  { 
    fTrackingAction
      ->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue)); 
  };   
}
