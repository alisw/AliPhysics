// $Id$
// Category: physics
//
// See the class description in the header file.

#include "AliEmptyPhysicsList.h"
#include "AliGlobals.h"

#include <G4Geantino.hh>
#include <G4ChargedGeantino.hh>

//_____________________________________________________________________________
AliEmptyPhysicsList::AliEmptyPhysicsList() {
//
  defaultCutValue = AliGlobals::DefaultCut();
  SetVerboseLevel(1);
}

//_____________________________________________________________________________
AliEmptyPhysicsList::~AliEmptyPhysicsList() {
//
}

// public methods

//_____________________________________________________________________________
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

//_____________________________________________________________________________
void AliEmptyPhysicsList::ConstructBosons()
{
// Constructs pseudo-particles only.
// ---

  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
}

//_____________________________________________________________________________
void AliEmptyPhysicsList::ConstructLeptons()
{
  // no leptons
}

//_____________________________________________________________________________
void AliEmptyPhysicsList::ConstructMesons()
{
 //  no mesons
}

//_____________________________________________________________________________
void AliEmptyPhysicsList::ConstructBarions()
{
 // no barions
}

//_____________________________________________________________________________
void AliEmptyPhysicsList::ConstructProcess()
{
// Constructs physics processes.
// ---

  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//_____________________________________________________________________________
void AliEmptyPhysicsList::ConstructEM()
{
  // no EM
}

//_____________________________________________________________________________
void AliEmptyPhysicsList::ConstructGeneral()
{
  // no Decay Process
}

//_____________________________________________________________________________
void AliEmptyPhysicsList::SetCuts()
{
// Sets the default range cut values for all defined particles.
// ---

  SetCutsWithDefault();
}

