
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// Some utility methods used for several wrapped functions
//
///////////////////////////////////////////////////////////////////

#ifndef WrapUtils_hh
#define WrapUtils_hh 1

//#include <iostream.h>
#include <iomanip>
#include "G4NavigationHistory.hh"
#include "G4ThreeVector.hh"

//Forward declarations

//StepAndLocation declaration. Used in:
// - WrapG1.cc
// - WrapNorml.cc
G4int StepAndLocation(G4ThreeVector &, const G4ThreeVector &, 
		      const G4double &, G4double &, G4double &, G4bool &,
		      G4bool &, const G4int &);

//EqualHistories declaration: true if histories are identical, otherwise false.
//Used in:
// - WrapG1.cc
bool EqualHistories(const G4NavigationHistory*, 
		    const G4NavigationHistory*);

// Common print utilities used in FGeometryInit.cc
inline G4std::ostream& setw10(G4std::ostream& os) { return os << G4std::setw(10);}
inline G4std::ostream& setscientific(G4std::ostream& os) { return os << G4std::setiosflags(G4std::ios::scientific);}
inline G4std::ostream& setfixed(G4std::ostream& os) { return os << G4std::setiosflags(G4std::ios::fixed);}
G4std::ostream& PrintHeader(G4std::ostream& os, const char* title);
G4String ToFlukaString(G4String str);

#endif
