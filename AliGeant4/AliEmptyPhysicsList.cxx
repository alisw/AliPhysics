// $Id$
// Category: physics
//
// See the class description in the header file.

#include "AliEmptyPhysicsList.h"
#include "AliGlobals.h"

#include <G4Geantino.hh>
#include <G4ChargedGeantino.hh>

AliEmptyPhysicsList::AliEmptyPhysicsList() {
//
  defaultCutValue = AliGlobals::DefaultCut();
  SetVerboseLevel(1);
}

AliEmptyPhysicsList::~AliEmptyPhysicsList() {
//
}

// public methods

void AliEmptyPhysicsList::ConstructParticle()
{
// In this method, static member functions should be called
// for all particles which you want to use.
// This ensures that objects of these particle types will be
// created in the program. 
// ---

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBarions();
}

void AliEmptyPhysicsList::ConstructBosons()
{
// Constructs pseudo-particles only.
// ---

  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
}

void AliEmptyPhysicsList::ConstructLeptons()
{
  // no leptons
}

void AliEmptyPhysicsList::ConstructMesons()
{
 //  no mesons
}

void AliEmptyPhysicsList::ConstructBarions()
{
 // no barions
}

void AliEmptyPhysicsList::ConstructProcess()
{
// Constructs physics processes.
// ---

  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

void AliEmptyPhysicsList::ConstructEM()
{
  // no EM
}

void AliEmptyPhysicsList::ConstructGeneral()
{
  // no Decay Process
}

void AliEmptyPhysicsList::SetCuts()
{
// Sets the default range cut values for all defined particles.
// ---

  SetCutsWithDefault();
}

