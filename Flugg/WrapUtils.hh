
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

// Commonly printed things used in FGeometryInit.cc
inline G4std::ostream& setw10(G4std::ostream& os) { return os << G4std::setw(10);}
inline G4std::ostream& setscientific(G4std::ostream& os) { return os << G4std::setiosflags(G4std::ios::scientific);}
inline G4std::ostream& setfixed(G4std::ostream& os) { return os << G4std::setiosflags(G4std::ios::fixed);}
G4std::ostream& PrintHeader(G4std::ostream& os, const char* title);
G4std::ostream& PrintMaterial(G4std::ostream& os, const char* title,
		       G4double Z, G4double A,
		       G4double density,
		       G4double index,
		       G4double N,
		       const char* name);
G4std::ostream& PrintCompound(G4std::ostream& os, const char* title,
		       G4int count,
		       const char* name,
		       G4double fraction,
		       G4double index);

#endif
