// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliEventActionMessenger.h"
#include "AliEventAction.h"
#include "AliGlobals.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithAnInteger.hh>

AliEventActionMessenger::AliEventActionMessenger(AliEventAction* eventAction)
  :fEventAction(eventAction)
{ 
//
  fEventDirectory = new G4UIdirectory("/aliEvent/");
  fEventDirectory->SetGuidance("AliEventAction control commands.");

  fDrawTracksCmd = new G4UIcmdWithAString("/aliEvent/drawTracks", this);
  fDrawTracksCmd->SetGuidance("Draw the tracks in the event");
  fDrawTracksCmd->SetGuidance("  Choice : NONE, CHARGED(default), ALL");
  fDrawTracksCmd->SetParameterName("Choice", true);
  fDrawTracksCmd->SetDefaultValue("CHARGED");
  fDrawTracksCmd->SetCandidates("NONE CHARGED ALL");
  fDrawTracksCmd->AvailableForStates(Idle);
 
  fVerboseCmd = new G4UIcmdWithAnInteger("/aliEvent/verbose", this);
  fVerboseCmd->SetGuidance("Set verbose level for AliEventAction");
  fVerboseCmd->SetParameterName("VerboseLevel", true);
  fVerboseCmd->SetDefaultValue(0);
  fVerboseCmd->SetRange("VerboseLevel >= 0 && VerboseLevel <= 2");
  fVerboseCmd->AvailableForStates(Idle);
}

AliEventActionMessenger::AliEventActionMessenger(){
//
}

AliEventActionMessenger::AliEventActionMessenger(
                                 const AliEventActionMessenger& right) {
//				 
  AliGlobals::Exception("AliEventActionMessenger is protected from copying.");
}

AliEventActionMessenger::~AliEventActionMessenger() {
//
  delete fEventDirectory;
  delete fDrawTracksCmd;
  delete fVerboseCmd;
}

// operators

AliEventActionMessenger& 
AliEventActionMessenger::operator=(const AliEventActionMessenger &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliEventActionMessenger is protected from assigning.");

  return *this;
}

// public methods

void AliEventActionMessenger::SetNewValue(G4UIcommand* command, 
       G4String newValue)
{ 
// Applies command to the associated object.
// ---

  if(command == fDrawTracksCmd)
  { 
    fEventAction->SetDrawFlag(newValue); 
  }   
  else if(command == fVerboseCmd)
  { 
    fEventAction
      ->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue)); 
  }   
}
