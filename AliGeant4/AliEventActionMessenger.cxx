// $Id$
// Category: event
//
// Author: I. Hrivnacova
//
// Class AliEventActionMessenger
// -----------------------------
// See the class description in the header file.

#include "AliEventActionMessenger.h"
#include "AliEventAction.h"
#include "AliGlobals.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithAnInteger.hh>

//_____________________________________________________________________________
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
}

//_____________________________________________________________________________
AliEventActionMessenger::AliEventActionMessenger(){
//
}

//_____________________________________________________________________________
AliEventActionMessenger::AliEventActionMessenger(
                                 const AliEventActionMessenger& right) {
//				 
  AliGlobals::Exception("AliEventActionMessenger is protected from copying.");
}

//_____________________________________________________________________________
AliEventActionMessenger::~AliEventActionMessenger() {
//
  delete fEventDirectory;
  delete fDrawTracksCmd;
}

// operators

//_____________________________________________________________________________
AliEventActionMessenger& 
AliEventActionMessenger::operator=(const AliEventActionMessenger &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliEventActionMessenger is protected from assigning.");

  return *this;
}

// public methods

//_____________________________________________________________________________
void AliEventActionMessenger::SetNewValue(G4UIcommand* command, 
       G4String newValue)
{ 
// Applies command to the associated object.
// ---

  if(command == fDrawTracksCmd)
  { 
    fEventAction->SetDrawFlag(newValue); 
  }   
}
