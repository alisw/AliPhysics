// $Id$
// Category: run
//
// See the class description in the header file.

#include "AliPrimaryGeneratorMessenger.h"
#include "AliPrimaryGeneratorAction.h"
#include "AliGlobals.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithAnInteger.hh>

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
  guidance = guidance + "  AliGenerator:      standard aliroot generator";
  fGeneratorCmd->SetGuidance(guidance);
  fGeneratorCmd->SetParameterName("Generator", true);
  fGeneratorCmd->SetCandidates("Gun Geantino ChargedGeantino AliGenerator");   
  fGeneratorCmd->SetDefaultValue("AliGenerator");
  fGeneratorCmd->AvailableForStates(PreInit,Idle);

  fNofParticlesCmd = new G4UIcmdWithAnInteger("/aliGenerator/nofParticles", this);
  fNofParticlesCmd->SetGuidance("Give number of primary particles:");
  fNofParticlesCmd->SetGuidance("It is applied only to \"gun\" type generators.");
  fNofParticlesCmd->SetParameterName("NofParticles", true);
  fNofParticlesCmd->SetDefaultValue(1);
  fNofParticlesCmd->SetRange("NofParticles >= 0");
  fNofParticlesCmd->AvailableForStates(PreInit,Idle);

  fVerboseCmd = new G4UIcmdWithAnInteger("/aliGenerator/verbose", this);
  fVerboseCmd->SetGuidance("Set verbose level for AliPrimaryGeneratorAction");
  fVerboseCmd->SetParameterName("VerboseLevel", true);
  fVerboseCmd->SetDefaultValue(0);
  fVerboseCmd->SetRange("VerboseLevel >= 0 && VerboseLevel <= 2");
  fVerboseCmd->AvailableForStates(Idle);
}

AliPrimaryGeneratorMessenger::AliPrimaryGeneratorMessenger() {
//
}

AliPrimaryGeneratorMessenger::AliPrimaryGeneratorMessenger(
                                 const AliPrimaryGeneratorMessenger& right) {
//				 
  AliGlobals::Exception(
    "AliPrimaryGeneratorMessenger is protected from copying.");
}

AliPrimaryGeneratorMessenger::~AliPrimaryGeneratorMessenger() {
//
  delete fPrimariesDirectory;
  delete fGeneratorCmd;
  delete fNofParticlesCmd;
  delete fVerboseCmd;
}

// operators

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
    else if (newValue == "AliGenerator")  
      fPrimaryGenAction->SetGenerator(kAliGenerator);       
  }
  else if( command == fNofParticlesCmd ) { 
    fPrimaryGenAction
      ->SetNofGunParticles(fNofParticlesCmd->GetNewIntValue(newValue)); 
  }   
  else if(command == fVerboseCmd) { 
    fPrimaryGenAction
      ->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue)); 
  }
}

