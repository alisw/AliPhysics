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

#ifndef G4ParticleDefinition_h
#define G4ParticleDefinition_h 1

#include "globals.hh"

class G4ParticleDefinition 
{
  public:
    G4ParticleDefinition();
    ~G4ParticleDefinition();

    const G4String& GetParticleName() const;
    const G4String& GetParticleType() const;
    G4double GetPDGCharge() const;

  private:
    G4String theParticleName;
    G4String theParticleType;
};

#endif

