// $Id$
// Category: physics
//
// Vector of control process flag values
// with convenient set/get methods

#ifndef TG4_FLAG_VECTOR_H
#define TG4_FLAG_VECTOR_H

#include "TG4Globals.h"
#include "TG3Flag.h"

#include <g4rw/tvordvec.h>

class G4VProcess;

class TG4FlagVector
{
  typedef G4RWTValOrderedVector<TG3FlagValue> TG3FlagVector;

  public:
    TG4FlagVector();
    TG4FlagVector(const TG4FlagVector& right);
    virtual ~TG4FlagVector();
    
    // operators
    TG4FlagVector& operator=(const TG4FlagVector& right);
    G4double operator[](G4int index) const;

    // set methods
    void SetG3Flag(TG3Flag g3Flag, G4double flagValue);
    void SetG3Defaults();
    
    // get methods
    G4int GetFlag(G4VProcess* process) const; 

  private:
    // data members
    TG3FlagVector*  fFlagVector; //vector of control process flag values
};

#endif //TG4_FLAG_VECTOR_H



