// $Id$
// Category: physics
//
// Constructor of electromagnetic physics.
// According to the corresponding part of:
// ExN04PhysicsList.hh, GEANT4 tag Name: geant4-01-01

#ifndef TG4_PHYSICS_CONSTRUCTOR_EM_H
#define TG4_PHYSICS_CONSTRUCTOR_EM_H

#include <G4VPhysicsConstructor.hh>
#include <globals.hh>

class TG4PhysicsConstructorEM: public G4VPhysicsConstructor
{
  public:
    TG4PhysicsConstructorEM(const G4String& name = "EM");
    virtual ~TG4PhysicsConstructorEM();

  protected:
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
};

#endif //TG4_PHYSICS_CONSTRUCTOR_H

