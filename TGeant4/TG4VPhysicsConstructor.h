// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorGeneral
// -----------------------------
// Base class for physics constructors with verbose.

#ifndef TG4_V_PHYSICS_CONSTRUCTOR_H
#define TG4_V_PHYSICS_CONSTRUCTOR_H

#include "TG4Verbose.h"

#include <G4VPhysicsConstructor.hh>
#include <globals.hh>

class TG4VPhysicsConstructor: public G4VPhysicsConstructor,
                              public TG4Verbose
{
  public:
    // TG4VPhysicsConstructor(); --> protected
    // TG4VPhysicsConstructor(const TG4VPhysicsConstructor& right);
    //  --> protected
    TG4VPhysicsConstructor(const G4String& name);
    TG4VPhysicsConstructor(const G4String& name, G4int verboseLevel);
    virtual ~TG4VPhysicsConstructor();

  protected:
    TG4VPhysicsConstructor();
    TG4VPhysicsConstructor(const TG4VPhysicsConstructor& right);

    // methods
          // construct particle and physics
    virtual void ConstructParticle() = 0;
    virtual void ConstructProcess() = 0;

    // overridden verbose methods
    virtual void  VerboseLevel(G4int level);
    virtual G4int VerboseLevel() const;
};

#endif //TG4_V_PHYSICS_CONSTRUCTOR_H

