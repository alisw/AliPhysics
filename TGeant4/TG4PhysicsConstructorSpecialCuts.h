// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorSpecialCuts
// --------------------------------------
// Constructor of special cuts.

#ifndef TG4_PHYSICS_CONSTRUCTOR_SPECIAL_CUTS_H
#define TG4_PHYSICS_CONSTRUCTOR_SPECIAL_CUTS_H

#include "TG4VPhysicsConstructor.h"

#include <globals.hh>

class TG4PhysicsConstructorSpecialCuts: public TG4VPhysicsConstructor
{
  public:
    TG4PhysicsConstructorSpecialCuts(const G4String& name = "Special cuts");
    TG4PhysicsConstructorSpecialCuts(G4int verboseLevel,
                                     const G4String& name = "Special cuts");
    virtual ~TG4PhysicsConstructorSpecialCuts();

  protected:
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
};

#endif //TG4_PHYSICS_CONSTRUCTOR_SPECIAL_CUTS_H

