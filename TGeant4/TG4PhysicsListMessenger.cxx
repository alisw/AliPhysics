// $Id$
// Category: physics
//
// See the class description in the header file.

#include "TG4PhysicsListMessenger.h"
#include "TG4PhysicsList.h"
#include "TG4Globals.h"

#include <G4UIcmdWithABool.hh>

TG4PhysicsListMessenger::TG4PhysicsListMessenger(TG4PhysicsList* physicsList)
  : fPhysicsList(physicsList)
{ 
//
  fSetOpticalCmd
     = new G4UIcmdWithABool("/g4mc/setOptical", this);
  fSetOpticalCmd->SetGuidance("Set Cerenkov and optical processes.");
  fSetOpticalCmd->SetParameterName("OpticalControl", false);
  fSetOpticalCmd->AvailableForStates(PreInit);

  fSetHadronCmd
     = new G4UIcmdWithABool("/g4mc/setHadron", this);
  fSetHadronCmd->SetGuidance("Set hadron processes.");
  fSetHadronCmd->SetParameterName("HadronControl", false);
  fSetHadronCmd->AvailableForStates(PreInit);

  fSetSpecialCutsCmd
     = new G4UIcmdWithABool("/g4mc/setSpecialCuts", this);
  fSetSpecialCutsCmd->SetGuidance("Set special cuts process.");
  fSetSpecialCutsCmd
    ->SetGuidance("!! Support for this option is under development.");
  fSetSpecialCutsCmd->SetParameterName("SpecialCutsControl", false);
  fSetSpecialCutsCmd->AvailableForStates(PreInit);

  fSetSpecialFlagsCmd
     = new G4UIcmdWithABool("/g4mc/setSpecialFlags", this);
  fSetSpecialFlagsCmd->SetGuidance("Set special flags process.");
  fSetSpecialFlagsCmd
    ->SetGuidance("!! Support for this option is under development.");
  fSetSpecialFlagsCmd->SetParameterName("SpecialFlagsControl", false);
  fSetSpecialFlagsCmd->AvailableForStates(PreInit);
}

TG4PhysicsListMessenger::TG4PhysicsListMessenger(
                            const TG4PhysicsListMessenger& right) {
//			     
  TG4Globals::Exception(
    "TG4PhysicsListMessenger is protected from copying.");
}

TG4PhysicsListMessenger::~TG4PhysicsListMessenger() {
//
  delete fSetOpticalCmd;
  delete fSetHadronCmd;
  delete fSetSpecialCutsCmd;
  delete fSetSpecialFlagsCmd;
}

// operators

TG4PhysicsListMessenger& 
TG4PhysicsListMessenger::operator=(const TG4PhysicsListMessenger& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "TG4PhysicsListMessenger is protected from assigning.");
    
  return *this;  
}    
          
// public methods

void TG4PhysicsListMessenger::SetNewValue(G4UIcommand* command, 
                                          G4String newValues)
{ 
// Applies command to the associated object.
// ---

  if (command == fSetOpticalCmd) {
    fPhysicsList
      ->SetOptical(fSetOpticalCmd->GetNewBoolValue(newValues)); 
  }    
  else if (command == fSetHadronCmd) {
    fPhysicsList
      ->SetHadron(fSetHadronCmd->GetNewBoolValue(newValues)); 
  }    
  else if (command == fSetSpecialCutsCmd) {
    fPhysicsList
      ->SetSpecialCuts(fSetSpecialCutsCmd->GetNewBoolValue(newValues)); 
  }    
  else if (command == fSetSpecialFlagsCmd) {
    fPhysicsList
      ->SetSpecialFlags(fSetSpecialFlagsCmd->GetNewBoolValue(newValues)); 
  }    
}
