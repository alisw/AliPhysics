// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4TG4PhysicsConstructorSpecialControls
// ---------------------------------------------
// Constructor of special controls of physics processes.

#ifndef TG4_PHYSICS_CONSTRUCTOR_SPECIAL_CONTROLS_H
#define TG4_PHYSICS_CONSTRUCTOR_SPECIAL_CONTROLS_H

#include "TG4VPhysicsConstructor.h"

#include <globals.hh>

class TG4PhysicsConstructorSpecialControls: public TG4VPhysicsConstructor
{
  public:
    TG4PhysicsConstructorSpecialControls(
      const G4String& name = "Special controls");
    TG4PhysicsConstructorSpecialControls(
      G4int verboseLevel,
      const G4String& name = "Special controls");
    virtual ~TG4PhysicsConstructorSpecialControls();

  protected:
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
};

#endif //TG4_PHYSICS_CONSTRUCTOR_SPECIAL_CONTROLS_H

