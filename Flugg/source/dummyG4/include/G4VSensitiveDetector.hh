// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id$
// GEANT4 tag $Name$
//

// dummy (required by G3toG4)

#ifndef G4VSensitiveDetector_h
#define G4VSensitiveDetector_h 1

#include "globals.hh"

class G4VSensitiveDetector 
{
  public:
    G4VSensitiveDetector();
    ~G4VSensitiveDetector();
    
    G4String GetName();
};




#endif

