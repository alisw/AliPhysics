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

// dummy (required by TGeant4)

#ifndef G4VProcess_h
#define G4VProcess_h 1

#include "globals.hh"

class G4VProcess 
{
  public:
    G4VProcess(const G4String&);
    ~G4VProcess();
    
    const G4String& GetProcessName() const {return theProcessName;}
      //  Returns the name of the process.

  private:
     G4String theProcessName;
};




#endif

