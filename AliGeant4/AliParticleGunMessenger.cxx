// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliParticleGunMessenger.h"
#include "AliParticleGun.h"
#include "AliGunParticle.h"
#include "AliGlobals.h"

#include <G4Geantino.hh>
#include <G4ThreeVector.hh>
#include <G4ParticleTable.hh>
#include <G4UIdirectory.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWith3Vector.hh>
#include <G4UIcmdWith3VectorAndUnit.hh>

AliParticleGunMessenger::AliParticleGunMessenger(AliParticleGun* gun)
  :fGun(gun)
{
//
  fParticleTable = G4ParticleTable::GetParticleTable();

  fGunDirectory = new G4UIdirectory("/aliGun/");
  fGunDirectory->SetGuidance("AliParticleGun control commands.");

  fListAvailableCmd 
    = new G4UIcmdWithoutParameter("/aliGun/listAvailable", this);
  fListAvailableCmd->SetGuidance("List available particles.");
  fListAvailableCmd->SetGuidance(" Invoke G4ParticleTable.");
  fListAvailableCmd->AvailableForStates(PreInit,Idle);

  fListCurrentCmd 
    = new G4UIcmdWithoutParameter("/aliGun/listCurrent", this);
  fListCurrentCmd->SetGuidance("List current particle properties.");
  fListCurrentCmd
    ->SetGuidance("(Use addParticle to add this particle to the gun.");
  fListCurrentCmd->AvailableForStates(PreInit,Idle);

  fParticleCmd 
    = new G4UIcmdWithAString("/aliGun/particle", this);
  fParticleCmd->SetGuidance("Set particle to be generated.");
  fParticleCmd->SetGuidance(" (geantino is default)");
  fParticleCmd->SetParameterName("particleName", true);
  fParticleCmd->SetDefaultValue("geantino");
  G4String candidateList; 
  G4int nofPTEntries = fParticleTable->entries();
  for (G4int i=0; i<nofPTEntries; i++)
  {
    candidateList += fParticleTable->GetParticleName(i);
    candidateList += " ";
  }
  fParticleCmd->SetCandidates(candidateList);
  fParticleCmd->AvailableForStates(PreInit,Idle);

  fMomentumCmd 
    = new G4UIcmdWith3VectorAndUnit("/aliGun/momentum", this);
  fMomentumCmd->SetGuidance("Set momentum.");
  fMomentumCmd->SetParameterName("Px","Py","Pz", true, true);
  fMomentumCmd->SetDefaultUnit("MeV");
  fMomentumCmd->SetUnitCategory("Energy"); 
  fMomentumCmd->SetRange("Px != 0 || Py != 0 || Pz != 0");
  fMomentumCmd->AvailableForStates(PreInit,Idle);

  fPositionCmd 
    = new G4UIcmdWith3VectorAndUnit("/aliGun/position", this);
  fPositionCmd->SetGuidance("Set starting position of the particle.");
  fPositionCmd->SetParameterName("X","Y","Z", true, true);
  fPositionCmd->SetDefaultUnit("cm");
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->AvailableForStates(PreInit,Idle);

  fTimeCmd 
    = new G4UIcmdWithADoubleAndUnit("/aliGun/time", this);
  fTimeCmd->SetGuidance("Set initial time of the particle.");
  fTimeCmd->SetParameterName("t0", true, true);
  fTimeCmd->SetDefaultUnit("ns");
  fTimeCmd->SetUnitCategory("Time");
  fTimeCmd->AvailableForStates(PreInit,Idle);
  
  fPolarizationCmd 
    = new G4UIcmdWith3Vector("/aliGun/polarization", this);
  fPolarizationCmd->SetGuidance("Set polarization.");
  fPolarizationCmd->SetParameterName("Px","Py","Pz", true, true); 
  fPolarizationCmd
    ->SetRange("Px>=-1. && Px<=1. && Py>=-1. && Py<=1. && Pz>=-1. && Pz<=1.");
  fPolarizationCmd->AvailableForStates(PreInit,Idle);
  
  fDirectionCmd 
    = new G4UIcmdWith3Vector("/aliGun/direction", this);
  fDirectionCmd->SetGuidance("Set momentum direction.");
  fDirectionCmd->SetGuidance("Direction needs not to be a unit vector.");
  fDirectionCmd->SetParameterName("Dx","Dy","Dz", true, true); 
  fDirectionCmd->SetRange("Dx != 0 || Dy != 0 || Dz != 0");
  fDirectionCmd->AvailableForStates(PreInit,Idle);
  
  fKinEnergyCmd 
    = new G4UIcmdWithADoubleAndUnit("/aliGun/kinEnergy", this);
  fKinEnergyCmd->SetGuidance("Set kinetic energy.");
  fKinEnergyCmd->SetParameterName("KineticEnergy", true, true);
  fKinEnergyCmd->SetDefaultUnit("GeV");
  fKinEnergyCmd->SetUnitCategory("Energy");
  fKinEnergyCmd->AvailableForStates(PreInit,Idle);

  fListCmd 
    = new G4UIcmdWithoutParameter("/aliGun/list",this);
  fListCmd->SetGuidance("List the Alice gun particles.");
  fListCmd->AvailableForStates(PreInit,Idle);

  fAddParticleCmd 
    = new G4UIcmdWithoutParameter("/aliGun/addParticle", this);
  fAddParticleCmd->SetGuidance("Add the particle to the Alice particle gun.");
  fAddParticleCmd->AvailableForStates(PreInit,Idle);

  fRemoveParticleCmd 
    = new G4UIcmdWithAnInteger("/aliGun/removeParticle", this);
  fRemoveParticleCmd
    ->SetGuidance("Remove the i-th particle friom the Alice particle gun.");
  fRemoveParticleCmd->SetParameterName("iParticle", false);
  fRemoveParticleCmd->SetRange("iParticle>=0");
  fRemoveParticleCmd->AvailableForStates(PreInit,Idle);

  fResetCmd 
    = new G4UIcmdWithoutParameter("/aliGun/reset", this);
  fResetCmd->SetGuidance("ReSet the Alice particle gun.");
  fResetCmd->AvailableForStates(PreInit,Idle);

  // Set initial value to AliGunParticle
  fParticle = new AliGunParticle();

  fParticle->SetParticleDefinition(G4Geantino::Geantino());
  fParticle->SetMomentumDirection(G4ThreeVector(1.0,0.0,0.0));
  fParticle->SetKineticEnergy(1.0*GeV);
  fParticle->SetPosition(G4ThreeVector(0.0*cm, 0.0*cm, 0.0*cm));
  fParticle->SetTime(0.0*ns);
  fParticle->SetPolarization(G4ThreeVector(0.0*cm, 0.0*cm, 0.0*cm));
}

AliParticleGunMessenger::AliParticleGunMessenger() {
//
}

AliParticleGunMessenger::AliParticleGunMessenger(
                                 const AliParticleGunMessenger& right) {
//				 
  AliGlobals::Exception("AliParticleGunMessenger is protected from copying.");
}

AliParticleGunMessenger::~AliParticleGunMessenger() {
//
  delete fListAvailableCmd;
  delete fParticleCmd;
  delete fMomentumCmd;
  delete fPositionCmd;
  delete fTimeCmd;
  delete fPolarizationCmd;
  delete fDirectionCmd;
  delete fKinEnergyCmd;
  delete fListCmd;
  delete fAddParticleCmd;
  delete fRemoveParticleCmd;
  delete fResetCmd;
  delete fGunDirectory;
  delete fParticle;
}

// operators

AliParticleGunMessenger& 
AliParticleGunMessenger::operator=(const AliParticleGunMessenger &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliParticleGunMessenger is protected from assigning.");

  return *this;
}

// public methods

void AliParticleGunMessenger::SetNewValue(G4UIcommand * command, 
       G4String newValues)
{
// Applies command to the associated object.
// ---

  // Alice gun particle set commands
  if (command==fListAvailableCmd)
  { 
    fParticleTable->DumpTable(); 
  } 
  else if (command == fParticleCmd)
  {
    G4ParticleDefinition* particleDef 
      = fParticleTable->FindParticle(newValues);
    if (particleDef != 0)
    { fParticle->SetParticleDefinition(particleDef); }
  } 
  else if (command == fListCurrentCmd)
  { 
    fParticle->Print(); 
  } 
  else if (command == fMomentumCmd)
  { 
    fParticle->SetMomentum(fMomentumCmd->GetNew3VectorValue(newValues)); 
  } 
  else if (command == fPositionCmd)
  { 
    fParticle->SetPosition(fDirectionCmd->GetNew3VectorValue(newValues)); 
  }
  else if (command == fTimeCmd)
  { 
    fParticle->SetTime(fTimeCmd->GetNewDoubleValue(newValues)); 
  }
  else if (command == fPolarizationCmd)
  { 
    fParticle
      ->SetPolarization(fPolarizationCmd->GetNew3VectorValue(newValues)); 
  } 
  else if (command == fDirectionCmd)
  { 
    fParticle
      ->SetMomentumDirection(fDirectionCmd->GetNew3VectorValue(newValues)); 
  } 
  else if (command == fKinEnergyCmd)
  { 
    fParticle->SetKineticEnergy(fKinEnergyCmd->GetNewDoubleValue(newValues)); 
  }

  // Alice particle gun commands 
  else if (command == fListCmd)
  { 
    fGun->List(); 
  }
  else if (command == fAddParticleCmd)
  { 
    fGun->AddParticle(fParticle);
    fParticle = new AliGunParticle(*fParticle); 
  } 
  else if (command == fRemoveParticleCmd)
  { 
    fGun->RemoveParticle(fRemoveParticleCmd->GetNewIntValue(newValues)); 
  }
  else if (command == fResetCmd)
  { 
    fGun->Reset(); 
  }
}

G4String AliParticleGunMessenger::GetCurrentValue(G4UIcommand * command)
{
// Returns current command parameters as string.
// ---

  G4String curValue;

  if( command==fDirectionCmd )
  { 
    curValue 
      = fDirectionCmd->ConvertToString(fParticle->GetMomentumDirection()); 
  } 
  else if( command==fKinEnergyCmd )
  { 
    curValue 
      = fKinEnergyCmd->ConvertToString(fParticle->GetKineticEnergy(),"GeV"); 
  } 
  else if( command==fPositionCmd )
  { 
    curValue = fPositionCmd->ConvertToString(fParticle->GetPosition(),"cm"); 
  } 
  else if( command==fTimeCmd )
  { 
    curValue = fTimeCmd->ConvertToString(fParticle->GetTime(),"ns"); 
  } 
  else if( command==fPolarizationCmd )
  { 
    curValue = fPolarizationCmd->ConvertToString(fParticle->GetPolarization()); 
  }
  
 return curValue;
}

