// $Id$
// Category: event
//
// See the class description in the header file.
 
#include "AliTrackingActionMessenger.h"
#include "AliTrackingAction.h"
#include "AliGlobals.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAnInteger.hh>

//_____________________________________________________________________________
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

  fNewVerboseCmd = new G4UIcmdWithAnInteger("/aliTracking/newVerbose", this);
  fNewVerboseCmd->SetGuidance("Set new verbose level (/tracking/verbose)");
  fNewVerboseCmd->SetGuidance("when a track with specified track ID ");
  fNewVerboseCmd->SetGuidance("(/aliTracking/newVerboseTrack)\n starts tracking");
  fNewVerboseCmd->SetParameterName("NewVerboseLevel", false);
  fNewVerboseCmd->SetRange("NewVerboseLevel >= 0 && NewVerboseLevel <= 5");
  fNewVerboseCmd->AvailableForStates(PreInit, Init, Idle);

  fNewVerboseTrackCmd = new G4UIcmdWithAnInteger("/aliTracking/newVerboseTrack", this);
  fNewVerboseTrackCmd->SetGuidance("Set the track ID for which the new verbose level");
  fNewVerboseTrackCmd->SetGuidance("(/aliTracking/newVerbose) will be applied.");
  fNewVerboseTrackCmd->SetParameterName("NewVerboseLevelTrackID", false);
  fNewVerboseTrackCmd->SetRange("NewVerboseLevelTrackID >= 0");
  fNewVerboseTrackCmd->AvailableForStates(PreInit, Init, Idle);
}

//_____________________________________________________________________________
AliTrackingActionMessenger::AliTrackingActionMessenger() {
//
}

//_____________________________________________________________________________
AliTrackingActionMessenger::AliTrackingActionMessenger(
                                 const AliTrackingActionMessenger& right) {
//				 
  AliGlobals::Exception(
    "AliTrackingActionMessenger is protected from copying.");
}

//_____________________________________________________________________________
AliTrackingActionMessenger::~AliTrackingActionMessenger() {
//
  delete fTrackingDirectory;
  delete fVerboseCmd;
  delete fNewVerboseCmd;
  delete fNewVerboseTrackCmd;
}

// operators

//_____________________________________________________________________________
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

//_____________________________________________________________________________
void AliTrackingActionMessenger::SetNewValue(G4UIcommand* command, 
       G4String newValue)
{ 
// Applies command to the associated object.
// ---

  if(command == fVerboseCmd) { 
    fTrackingAction
      ->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue)); 
  }   
  else if(command == fNewVerboseCmd) { 
    fTrackingAction
      ->SetNewVerboseLevel(fNewVerboseCmd->GetNewIntValue(newValue)); 
  }   
  else if(command == fNewVerboseTrackCmd) { 
    fTrackingAction
      ->SetNewVerboseTrackID(fNewVerboseTrackCmd->GetNewIntValue(newValue)); 
  }   
}
