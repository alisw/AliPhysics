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

// dummy (required by TGeant4)

#ifndef G4Track_h
#define G4Track_h 1

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"

class G4Track 
{
  public:
    G4Track();
    ~G4Track();

    G4ParticleDefinition* GetDefinition() const {return 0;}
    const G4VProcess* GetCreatorProcess() const {return 0;}
};

#endif

