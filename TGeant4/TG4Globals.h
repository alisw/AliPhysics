// $Id$
// Category: global
//
// Class for generally used basic types and functions.
// It is protected from instantiating (only static data members
// and static methods are defined).

#ifndef TG4_GLOBALS_H
#define TG4_GLOBALS_H

#include <globals.hh>
#include <g4std/vector>
#include <g4std/set>
#include <g4rw/tvordvec.h>
#include <g4rw/tpordvec.h>
#include <G4RotationMatrix.hh>

class G4Material;
class G4Element;

typedef G4std::vector<G4bool>  TG4boolVector;
typedef G4std::vector<G4int>   TG4intVector;
typedef G4std::vector<const G4RotationMatrix*> TG4RotationMatrixVector;
typedef G4std::set <G4String, G4std::less<G4String> > TG4StringSet; 
typedef G4RWTValOrderedVector<G4double> TG4doubleVector;
typedef G4RWTValOrderedVector<G4String> TG4StringVector;
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
    static void AppendNumberToString(G4String& string, G4int number);

  protected:
    TG4Globals();  
      // only typedefs's and static methods
};  

#endif //ALGLOBALS_H
