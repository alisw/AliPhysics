// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliSteppingActionMessenger.h"
#include "AliSteppingAction.h"
#include "AliGlobals.h"

#include <G4UIcmdWithAnInteger.hh>

AliSteppingActionMessenger::AliSteppingActionMessenger(
                              AliSteppingAction* trackingAction)
  :fSteppingAction(trackingAction)
{
// 
  fLoopVerboseCmd = new G4UIcmdWithAnInteger("/aliTracking/loopVerbose", this);
  fLoopVerboseCmd
    ->SetGuidance("Set tracking verbose level for detected looping tracks.");
  fLoopVerboseCmd->SetParameterName("LoopVerboseLevel", true);
  fLoopVerboseCmd->SetDefaultValue(1);
  fLoopVerboseCmd->SetRange("LoopVerboseLevel >= 0 && LoopVerboseLevel <= 5");
  fLoopVerboseCmd->AvailableForStates(Idle);

  fMaxNofStepsCmd = new G4UIcmdWithAnInteger("/aliTracking/maxNofSteps", this);
  fMaxNofStepsCmd
    ->SetGuidance("Set tracking verbose level for detected looping tracks.");
  fMaxNofStepsCmd->SetParameterName("MaxNofSteps", false);
  fMaxNofStepsCmd->SetRange("MaxNofSteps >= 0");
  fMaxNofStepsCmd->AvailableForStates(Idle);
}

AliSteppingActionMessenger::AliSteppingActionMessenger() {
//
}

AliSteppingActionMessenger::AliSteppingActionMessenger(
                                 const AliSteppingActionMessenger& right) {
//				 
  AliGlobals::Exception(
    "AliSteppingActionMessenger is protected from copying.");
}

AliSteppingActionMessenger::~AliSteppingActionMessenger() {
//
  delete fLoopVerboseCmd;
  delete fMaxNofStepsCmd;
}

// operators

AliSteppingActionMessenger& 
AliSteppingActionMessenger::operator=(const AliSteppingActionMessenger &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception(
    "AliSteppingActionMessenger is protected from assigning.");

  return *this;
}

// public methods

void AliSteppingActionMessenger::SetNewValue(G4UIcommand* command, 
       G4String newValue)
{ 
// Applies command to the associated object.
// ---

  if(command == fLoopVerboseCmd) { 
    fSteppingAction
      ->SetLoopVerboseLevel(fLoopVerboseCmd->GetNewIntValue(newValue)); 
  }   
  else if(command == fMaxNofStepsCmd) { 
    fSteppingAction
      ->SetMaxNofSteps(fMaxNofStepsCmd->GetNewIntValue(newValue)); 
  }   
}
