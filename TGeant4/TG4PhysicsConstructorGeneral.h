// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorGeneral
// -----------------------------
// Constructor of physics for ions.
// According to ExN04IonPhysics.hh, GEANT4 tag Name: geant4-03-02

#ifndef TG4_PHYSICS_CONSTRUCTOR_GENERAL_H
#define TG4_PHYSICS_CONSTRUCTOR_GENERAL_H

#include <G4VPhysicsConstructor.hh>
#include <G4Decay.hh>
#include <globals.hh>

class TG4PhysicsConstructorGeneral: public G4VPhysicsConstructor
{
  public:
    TG4PhysicsConstructorGeneral(const G4String& name = "General");
    virtual ~TG4PhysicsConstructorGeneral();

  protected:
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    // data members
    G4Decay fDecayProcess; // decay process
};

#endif //TG4_PHYSICS_CONSTRUCTOR_GENERAL_H

