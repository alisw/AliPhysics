// $Id$
// Category: run
//
// Author: I. Hrivnacova
//
// Class AliPrimaryGeneratorMessenger
// ----------------------------------
// See the class description in the header file.

#include "AliPrimaryGeneratorMessenger.h"
#include "AliPrimaryGeneratorAction.h"
#include "AliGlobals.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithAnInteger.hh>

//_____________________________________________________________________________
AliPrimaryGeneratorMessenger::AliPrimaryGeneratorMessenger(
  AliPrimaryGeneratorAction* primaryGenAction)
  : fPrimaryGenAction(primaryGenAction)
{
// 
  fPrimariesDirectory = new G4UIdirectory("/aliGenerator/");
  fPrimariesDirectory->SetGuidance("AliPrimaryGeneratorAction control commands.");

  fGeneratorCmd = new G4UIcmdWithAString("/aliGenerator/set", this);
  G4String guidance =   "Set one of the defined primary generators:\n";
  guidance = guidance + "  Gun:               particle gun (default)\n";
  guidance = guidance + "  Geantino:          geantino with random momentum(default)\n";
  guidance = guidance + "  ChargedGeantino:   chargedgeantino with random momentum\n";
  guidance = guidance + "  Stack:             standard generator from MC stack";
  fGeneratorCmd->SetGuidance(guidance);
  fGeneratorCmd->SetParameterName("Generator", true);
  fGeneratorCmd->SetCandidates("Gun Geantino ChargedGeantino Stack");   
  fGeneratorCmd->SetDefaultValue("Stack");
  fGeneratorCmd->AvailableForStates(PreInit,Idle);

  fNofParticlesCmd = new G4UIcmdWithAnInteger("/aliGenerator/nofParticles", this);
  fNofParticlesCmd->SetGuidance("Give number of primary particles:");
  fNofParticlesCmd->SetGuidance("It is applied only to \"gun\" type generators.");
  fNofParticlesCmd->SetParameterName("NofParticles", true);
  fNofParticlesCmd->SetDefaultValue(1);
  fNofParticlesCmd->SetRange("NofParticles >= 0");
  fNofParticlesCmd->AvailableForStates(PreInit,Idle);
}

//_____________________________________________________________________________
AliPrimaryGeneratorMessenger::AliPrimaryGeneratorMessenger() {
//
}

//_____________________________________________________________________________
AliPrimaryGeneratorMessenger::AliPrimaryGeneratorMessenger(
                                 const AliPrimaryGeneratorMessenger& right) {
//				 
  AliGlobals::Exception(
    "AliPrimaryGeneratorMessenger is protected from copying.");
}

//_____________________________________________________________________________
AliPrimaryGeneratorMessenger::~AliPrimaryGeneratorMessenger() {
//
  delete fPrimariesDirectory;
  delete fGeneratorCmd;
  delete fNofParticlesCmd;
}

// operators

//_____________________________________________________________________________
AliPrimaryGeneratorMessenger& 
AliPrimaryGeneratorMessenger::operator=(
                                const AliPrimaryGeneratorMessenger &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception(
    "AliPrimaryGeneratorMessenger is protected from assigning.");

  return *this;
}

// public methods

//_____________________________________________________________________________
void AliPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,
       G4String newValue)
{ 
// Applies command to the associated object.
// ---

  if( command == fGeneratorCmd ) { 
    if (newValue == "Gun") 
      fPrimaryGenAction->SetGenerator(kGun); 
    else if (newValue == "Geantino") 
      fPrimaryGenAction->SetGenerator(kGeantino); 
    else if (newValue == "ChargedGeantino")  
      fPrimaryGenAction->SetGenerator(kChargedGeantino); 
    else if (newValue == "Stack")  
      fPrimaryGenAction->SetGenerator(kStack);       
  }
  else if( command == fNofParticlesCmd ) { 
    fPrimaryGenAction
      ->SetNofGunParticles(fNofParticlesCmd->GetNewIntValue(newValue)); 
  }   
}

