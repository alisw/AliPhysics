// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4Globals
// ----------------
// Class provides the basic types and functions of general use. 
// It is protected from instantiating (only static data members
// and static methods are defined).

#ifndef TG4_GLOBALS_H
#define TG4_GLOBALS_H

#include "TG4G3Control.h"

#include <globals.hh>
#include <g4std/vector>
#include <g4std/set>
#include <G4RotationMatrix.hh>

class G4Material;
class G4Element;

// basic types containers
typedef G4std::vector<G4bool>   TG4boolVector;
typedef G4std::vector<G4int>    TG4intVector;
typedef G4std::vector<G4double> TG4doubleVector;
typedef G4std::vector<G4String> TG4StringVector;
typedef G4std::set <G4String, G4std::less<G4String> > TG4StringSet; 

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
    static G4bool Compare(G4bool activation, TG4G3ControlValue controlValue);
    static void PrintStars(G4bool emptyLineFirst);

  protected:
    TG4Globals();  
      // only typedefs's and static methods
};  

#endif //ALGLOBALS_H
