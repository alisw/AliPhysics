
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// Some utility methods used for several wrapped functions
//
///////////////////////////////////////////////////////////////////

#ifndef WrapUtils_hh
#define WrapUtils_hh 1

#include "G4NavigationHistory.hh"
#include "G4ThreeVector.hh"
#include <iostream.h>
#include <iomanip.h>

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
inline ostream& setw10(ostream& os) { return os << G4std::setw(10);}
inline ostream& setscientific(ostream& os) { return os << G4std::setiosflags(G4std::ios::scientific);}
inline ostream& setfixed(ostream& os) { return os << G4std::setiosflags(G4std::ios::fixed);}
ostream& PrintHeader(ostream& os, const char* title);
ostream& PrintMaterial(ostream& os, const char* title,
		       G4double Z, G4double A,
		       G4double density,
		       G4double index,
		       G4double N,
		       const char* name);
ostream& PrintCompound(ostream& os, const char* title,
		       G4int count,
		       const char* name,
		       G4double fraction,
		       G4double index);

#endif
