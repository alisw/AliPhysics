// $Id$
// Category: physics
//
// Constructor of special controls of physics processes.

#ifndef TG4_PHYSICS_CONSTRUCTOR_SPECIAL_CONTROLS_H
#define TG4_PHYSICS_CONSTRUCTOR_SPECIAL_CONTROLS_H

#include <G4VPhysicsConstructor.hh>
#include <globals.hh>

class TG4PhysicsConstructorSpecialControls: public G4VPhysicsConstructor
{
  public:
    TG4PhysicsConstructorSpecialControls(
      const G4String& name = "Special controls");
    virtual ~TG4PhysicsConstructorSpecialControls();

  protected:
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
};

#endif //TG4_PHYSICS_CONSTRUCTOR_SPECIAL_CONTROLS_H

