// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4VPhysicsConstructor
// ------------------------------
// See the class description in the header file.
// According to ExN04IonPhysics.cc,v 1.1.2.1 2001/06/28 19:07:37 gunter Exp 
// GEANT4 tag Name: geant4-03-02

#include "TG4VPhysicsConstructor.h"
#include "TG4Globals.h"

//_____________________________________________________________________________
TG4VPhysicsConstructor::TG4VPhysicsConstructor(const G4String& name)
  : G4VPhysicsConstructor(name),
    TG4Verbose(G4String("physics" + name))
{
//
  VerboseLevel(1);
}

//_____________________________________________________________________________
TG4VPhysicsConstructor::TG4VPhysicsConstructor(const G4String& name,
                                               G4int verboseLevel)
  : G4VPhysicsConstructor(name),
    TG4Verbose(G4String("physics" + name))
{
//
  VerboseLevel(verboseLevel);
}

//_____________________________________________________________________________
TG4VPhysicsConstructor::TG4VPhysicsConstructor(
                                     const TG4VPhysicsConstructor& right)    
  : TG4Verbose("") {
//
  TG4Globals::Exception("TG4VPhysicsConstructor is protected from copying.");
}

//_____________________________________________________________________________
TG4VPhysicsConstructor::TG4VPhysicsConstructor()  
  : TG4Verbose("") {
//
}

//_____________________________________________________________________________
TG4VPhysicsConstructor::~TG4VPhysicsConstructor() {
//
}

//
// public methods
//

//_____________________________________________________________________________
void TG4VPhysicsConstructor::VerboseLevel(G4int level)
{
// Sets the same value to G4VPhysicsConstructor verbose
// level and TG4Verbose level.
// ---

   TG4Verbose::VerboseLevel(level);
   
   // verbose in G4VPhysicsConstructor 
   SetVerboseLevel(level);
}


//_____________________________________________________________________________
G4int TG4VPhysicsConstructor::VerboseLevel() const
{
// Returns TG4Verbose level value.
// ---

   return TG4Verbose::VerboseLevel();
}
