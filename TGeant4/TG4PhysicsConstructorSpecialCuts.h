// $Id$
// Category: physics
//
// Constructor of special cuts.

#ifndef TG4_PHYSICS_CONSTRUCTOR_SPECIAL_CUTS_H
#define TG4_PHYSICS_CONSTRUCTOR_SPECIAL_CUTS_H

#include <G4VPhysicsConstructor.hh>
#include <globals.hh>

class TG4PhysicsConstructorSpecialCuts: public G4VPhysicsConstructor
{
  public:
    TG4PhysicsConstructorSpecialCuts(const G4String& name = "Special cuts");
    virtual ~TG4PhysicsConstructorSpecialCuts();

  protected:
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
};

#endif //TG4_PHYSICS_CONSTRUCTOR_SPECIAL_CUTS_H

