
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// WrapDN.hh - Sara Vanini - 17/XI/98
//
// Wrapper for setting DNEAR option on fluka side. Must return 0 
// if user doesn't want Fluka to use DNEAR to compute the 
// step (the same effect is obtained with the GLOBAL (WHAT(3)=-1)
// card in fluka input), returns 1 if user wants Fluka always to 
// use DNEAR (in this case, be sure that GEANT4 DNEAR is unique, 
// coming from all directions!!!).
//
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
///////////////////////////////////////////////////////////////////

//#ifndef idnrwr
//#define idnrwr idnrwr_

#include "Wrappers.hh"
#include "globals.hh"

G4int idnrwr(const G4int & nreg, const G4int & mlat) 

{
//flag
#ifdef G4GEOMETRY_DEBUG
	G4cout<<"================== IDNRWR ================="<<G4endl;
#endif 

// returns 0 if user doesn't want Fluka to use DNEAR to compute the 
// step (the same effect is obtained with the GLOBAL (WHAT(3)=-1)
// card in fluka input), returns 1 if (be sure that GEANT4 DNEAR is unique, 
// coming from all directions!!!) user wants Fluka always to use DNEAR.

return 0;
}
//#endif
