// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliDetSwitchVectorMessenger
// ------------------------------------
// See the class description in the header file.

#include "AliDetSwitchVectorMessenger.h"
#include "AliDetSwitchVector.h"
#include "AliGlobals.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithoutParameter.hh>

//_____________________________________________________________________________
AliDetSwitchVectorMessenger::AliDetSwitchVectorMessenger(
                                         AliDetSwitchVector* switchVector)
  : fDetSwitchVector(switchVector)
{
//
  fSwitchOnCmd = new G4UIcmdWithAString("/aliDet/switchOn", this);
  fSwitchOnCmd->SetGuidance("Define the module to be built.");
  fSwitchOnCmd->SetGuidance("Available modules:");
  G4String listAvailableDets = "NONE, ALL";  
  fSwitchOnCmd->SetGuidance(listAvailableDets);
  fSwitchOnCmd->SetParameterName("module", false);
  fSwitchOnCmd->AvailableForStates(PreInit);;

  fSwitchOffCmd = new G4UIcmdWithAString("/aliDet/switchOff", this);
  fSwitchOffCmd->SetGuidance("Define the module not to be built.");
  fSwitchOffCmd->SetGuidance("Available modules:");
  G4String listDetsNames = "ALL"; 
  fSwitchOffCmd->SetGuidance(listDetsNames);
  fSwitchOffCmd->SetParameterName("module", false);
  fSwitchOffCmd->AvailableForStates(PreInit);;

  fListCmd 
    = new G4UIcmdWithoutParameter("/aliDet/list", this);
  fListCmd->SetGuidance("List the currently switched modules.");
  fListCmd
    ->AvailableForStates(PreInit, Init, Idle, GeomClosed, EventProc);

  fListAvailableCmd 
    = new G4UIcmdWithoutParameter("/aliDet/listAvailable", this);
  fListAvailableCmd->SetGuidance("List all available modules.");
  fListAvailableCmd
    ->AvailableForStates(PreInit, Init, Idle, GeomClosed, EventProc);
}

//_____________________________________________________________________________
AliDetSwitchVectorMessenger::AliDetSwitchVectorMessenger() {
//
}

//_____________________________________________________________________________
AliDetSwitchVectorMessenger::AliDetSwitchVectorMessenger(
                                   const AliDetSwitchVectorMessenger& right)
{
//
  AliGlobals::Exception(
    "AliDetSwitchVectorMessenger is protected from copying.");
}

//_____________________________________________________________________________
AliDetSwitchVectorMessenger::~AliDetSwitchVectorMessenger() {
//
  delete fSwitchOnCmd;
  delete fSwitchOffCmd;
  delete fListCmd;
  delete fListAvailableCmd;
}

// operators

//_____________________________________________________________________________
AliDetSwitchVectorMessenger& 
AliDetSwitchVectorMessenger::operator=(const AliDetSwitchVectorMessenger& right)
{
  // check assignement to self
  if (this == &right) return *this;

  AliGlobals::Exception(
     "AliDetSwitchVectorMessenger is protected from assigning.");
    
  return *this;  
}    
          
// private methods

//_____________________________________________________________________________
void AliDetSwitchVectorMessenger::SetGuidance() 
{
// Updates guidance, candidates list.
// ---

  G4String listAvailableDets = "NONE, ALL, ";  
  listAvailableDets 
    = listAvailableDets + fDetSwitchVector->GetAvailableDetsListWithCommas();
  fSwitchOnCmd->SetGuidance(listAvailableDets);

  G4String listDetsNames = "ALL, "; 
  listDetsNames
    = listDetsNames + fDetSwitchVector->GetDetNamesListWithCommas();
  fSwitchOffCmd->SetGuidance(listDetsNames);
}

//_____________________________________________________________________________
void AliDetSwitchVectorMessenger::SetCandidates() 
{
// Builds candidates list.
// ---

  G4String candidatesList = "NONE ALL ";
  candidatesList += fDetSwitchVector->GetDetNamesList();;
  candidatesList += fDetSwitchVector->GetAvailableDetsList();
  fSwitchOnCmd->SetCandidates(candidatesList);

  candidatesList = "ALL ";
  candidatesList += fDetSwitchVector->GetDetNamesList();;
  fSwitchOffCmd->SetCandidates(candidatesList);
}

// public methods
  
//_____________________________________________________________________________
void AliDetSwitchVectorMessenger::SetNewValue(G4UIcommand* command, 
                                              G4String newValues)
{
// Applies command to the associated object.
// ---

  if (command == fSwitchOnCmd) {  
    fDetSwitchVector->SwitchDetOn(newValues); 
  }
  else if (command == fSwitchOffCmd) {  
    fDetSwitchVector->SwitchDetOff(newValues); 
  }
  else if (command == fListCmd) {  
    fDetSwitchVector->PrintSwitchedDets(); 
  }
  else if (command == fListAvailableCmd) {  
    fDetSwitchVector->PrintAvailableDets(); 
  }
}

//_____________________________________________________________________________
void AliDetSwitchVectorMessenger::Update() 
{
// Updates guidance, candidates list.
// ---

  SetGuidance();
  SetCandidates();

  // set default values to a detector
  fDetSwitchVector->SwitchDetOn("NONE");
}
