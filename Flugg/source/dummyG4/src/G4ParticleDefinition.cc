// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id$
// GEANT4 tag $Name$
//

// dummy (required by G3toG4)

#include "G4ParticleDefinition.hh"

G4ParticleDefinition::G4ParticleDefinition()
  : theParticleName(""),
    theParticleType("")
{
}

G4ParticleDefinition::~G4ParticleDefinition()
{
}

const G4String& G4ParticleDefinition::GetParticleName() const
{ 
  return theParticleName;
}

const G4String& G4ParticleDefinition::GetParticleType() const 
{ 
  return theParticleType;
}

G4double G4ParticleDefinition::GetPDGCharge() const 
{ 
  return 0; 
}
