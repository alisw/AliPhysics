// $Id$
// Category: physics
//
// Constructor of optical physics.
// According to ExN06PhysicsList (geant4 1.1)

#ifndef TG4_PHYSICS_CONSTRUCTOR_OPTICAL_H
#define TG4_PHYSICS_CONSTRUCTOR_OPTICAL_H

#include <G4VPhysicsConstructor.hh>
#include <globals.hh>

class TG4PhysicsConstructorOptical: public G4VPhysicsConstructor
{
  public:
    TG4PhysicsConstructorOptical(const G4String& name = "Optical");
    virtual ~TG4PhysicsConstructorOptical();

  protected:
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
};

#endif //TG4_PHYSICS_CONSTRUCTOR_OPTICAL_H

