// $Id$
// Category: physics
//
// Constructor of hadron physics.
// According to corresponding part of:
// ExN04PhysicsList.hh, GEANT4 tag Name: geant4-01-01

#ifndef TG4_PHYSICS_CONSTRUCTOR_HADRON_H
#define TG4_PHYSICS_CONSTRUCTOR_HADRON_H

#include <G4VPhysicsConstructor.hh>
#include <globals.hh>

class TG4PhysicsConstructorHadron: public G4VPhysicsConstructor
{
  public:
    TG4PhysicsConstructorHadron(const G4String& name = "Hadron");
    virtual ~TG4PhysicsConstructorHadron();

  protected:
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
};

#endif //TG4_PHYSICS_CONSTRUCTOR_HADRON_H

