
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// WrapG1.hh - Sara Vanini
//
// Dummy wrapper (in fluka: for geometry debugging)
//
//////////////////////////////////////////////////////////////////

//#ifndef lkdbwr
//#define lkdbwr lkdbwr_

#include "Wrappers.hh"
#include "globals.hh"

void lkdbwr(G4double& pSx, G4double& pSy, G4double& pSz,
	    G4double* pV, const G4int& oldReg, const G4int& oldLttc,
	    G4int& newReg, G4int& flagErr, G4int& newLttc)
{
  //flag
#ifdef G4GEOMETRY_DEBUG
  G4cout<<"============= LKDBWR =============="<<G4endl;
#endif
  
  //return region number and dummy variables
  newReg=0;   
  newLttc=0;
  flagErr=-1; 
}
//#endif
