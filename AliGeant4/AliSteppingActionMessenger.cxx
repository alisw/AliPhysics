// $Id$
// Category: event
//
// Author: I. Hrivnacova
//
// Class AliSteppingActionMessenger
// --------------------------------
// See the class description in the header file.

#include "AliSteppingActionMessenger.h"
#include "AliSteppingAction.h"
#include "AliGlobals.h"

#include <G4UIcmdWithAnInteger.hh>

//_____________________________________________________________________________
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
    ->SetGuidance("Set maximum number of steps allowed.");
  fMaxNofStepsCmd->SetParameterName("MaxNofSteps", false);
  fMaxNofStepsCmd->SetRange("MaxNofSteps >= 0");
  fMaxNofStepsCmd->AvailableForStates(Idle);
}

//_____________________________________________________________________________
AliSteppingActionMessenger::AliSteppingActionMessenger() {
//
}

//_____________________________________________________________________________
AliSteppingActionMessenger::AliSteppingActionMessenger(
                                 const AliSteppingActionMessenger& right) {
//				 
  AliGlobals::Exception(
    "AliSteppingActionMessenger is protected from copying.");
}

//_____________________________________________________________________________
AliSteppingActionMessenger::~AliSteppingActionMessenger() {
//
  delete fLoopVerboseCmd;
  delete fMaxNofStepsCmd;
}

// operators

//_____________________________________________________________________________
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

//_____________________________________________________________________________
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
