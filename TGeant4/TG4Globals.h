// $Id$
// Category: global
//
// Class for generally used basic types and functions.
// It is protected from instantiating (only static data members
// and static methods are defined).

#ifndef TG4_GLOBALS_H
#define TG4_GLOBALS_H

#include "TG3Flag.h"

#include <globals.hh>

#include <g4rw/tvordvec.h>
#include <g4rw/tpordvec.h>

#include <vector.h>

class G4Material;
class G4Element;

//typedef G4RWTValOrderedVector<G4bool>   TG4boolVector;
typedef vector<G4bool>   TG4boolVector;
typedef G4RWTValOrderedVector<G4double> TG4doubleVector;
typedef G4RWTValOrderedVector<G4String> TG4StringVector;
typedef G4RWTValOrderedVector<TG3FlagValue> TG3FlagVector;
typedef G4RWTPtrOrderedVector<G4Material>   TG4MaterialVector;
typedef G4RWTPtrOrderedVector<G4Element>    TG4ElementVector;  

class TG4Globals
{
  public:
    // --> protected 
    // TG4Globals();
    virtual ~TG4Globals();

    // static methods
    static void Exception(const char* string=0);
      // Global error function prints string to cerr, and aborts
      // program - according to G4Exception.cc
    static void Warning(const char* string=0);
      // Global warning function prints string to cerr

  protected:
    TG4Globals();  
      // only typedefs's and static methods
};  

#endif //ALGLOBALS_H
