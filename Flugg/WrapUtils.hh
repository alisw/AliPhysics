
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
inline std::ostream& setw10(std::ostream& os) { return os << std::setw(10);}
inline std::ostream& setscientific(std::ostream& os) { return os << std::setiosflags(std::ios::scientific);}
inline std::ostream& setfixed(std::ostream& os) { return os << std::setiosflags(std::ios::fixed);}
std::ostream& PrintHeader(std::ostream& os, const char* title);
G4String ToFlukaString(G4String str);

#endif
