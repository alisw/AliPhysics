
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// Interface for Flugg Wrappers
//
///////////////////////////////////////////////////////////////////

#ifndef WRAPPERS_HH
#define WRAPPERS_HH

#include "globals.hh"

#define idnrwr idnrwr_
#define g1wr g1wr_
#define g1rtwr g1rtwr_
#define conhwr conhwr_
#define inihwr inihwr_
#define jomiwr jomiwr_
#define lkdbwr lkdbwr_
#define lkfxwr lkfxwr_
#define lkmgwr lkmgwr_
#define lkwr lkwr_
#define magfld magfld_
#define nrmlwr nrmlwr_
#define rgrpwr rgrpwr_
#define isvhwr isvhwr_


#define G4GEOMETRY_DEBUG 1

// WrapDN

extern "C" G4int idnrwr(const G4int & nreg, const G4int & mlat);

// WrapG1

extern "C" void  g1wr(G4double& pSx, G4double& pSy, G4double& pSz, G4double* pV,
                      G4int& oldReg, const G4int& oldLttc, G4double& propStep,
                      G4int& nascFlag, G4double& retStep, G4int& newReg,
	              G4double& saf, G4int& newLttc, G4int& LttcFlag,
                      G4double* sLt, G4int* jrLt);

// WrapG1RT

extern "C" void g1rtwr(void);

// WrapIncrHist

extern "C" void conhwr(G4int& intHist, G4int* incrCount); 

// WrapIniHist

extern "C" void inihwr(G4int& intHist);                   

// WrapInit

extern "C" void jomiwr(const G4int & nge, const G4int& lin, const G4int& lou,
                       G4int& flukaReg);

// WrapLookDB

extern "C" void lkdbwr(G4double& pSx, G4double& pSy, G4double& pSz,
                       G4double* pV, const G4int& oldReg, const G4int& oldLttc,
               	       G4int& newReg, G4int& flagErr, G4int& newLttc);

// WrapLookFX

extern "C" void lkfxwr(G4double& pSx, G4double& pSy, G4double& pSz,
                       G4double* pV, const G4int& oldReg, const G4int& oldLttc,
                       G4int& newReg, G4int& flagErr, G4int& newLttc);
	    
// WrapLookMG

extern "C" void lkmgwr(G4double& pSx, G4double& pSy, G4double& pSz,
                       G4double* pV, const G4int& oldReg, const G4int& oldLttc,
		       G4int& flagErr, G4int& newReg, G4int& newLttc);
	    
// WrapLookZ

extern "C" void lkwr(G4double& pSx, G4double& pSy, G4double& pSz,
                     G4double* pV, const G4int& oldReg, const G4int& oldLttc,
	             G4int& newReg, G4int& flagErr, G4int& newLttc);

// WrapMag

extern "C" void magfld(const G4double& pX, const G4double& pY, const G4double& pZ,
                       G4double& cosBx, G4double& cosBy, G4double& cosBz, 
                       G4double& Bmag, G4int& reg, G4int& idiscflag);
	    
// WrapNorml

extern "C" void nrmlwr(G4double& pSx, G4double& pSy, G4double& pSz,
                       G4double& pVx, G4double& pVy, G4double& pVz,
	               G4double* norml, const G4int& oldReg, 
	               const G4int& newReg, G4int& flagErr);

// WrapReg

extern "C" void rgrpwr(const G4int& flukaReg, const G4int& ptrLttc, G4int& g4Reg,
                       G4int* indMother, G4int* repMother, G4int& depthFluka);

// WrapSavHist
	    
extern "C" G4int isvhwr(const G4int& fCheck, const G4int& intHist);

#endif //WRAPPERS_HH

