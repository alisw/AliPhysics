
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// WrapLookFX.hh - Sara Vanini - 24/III/00
//
// Wrapper for localisation of particle to fix particular conditions.
// At the moment is the same as WrapLookZ.hh. 
//
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
// modified 17/06/02: by I. Gonzalez. STL migration
//
//////////////////////////////////////////////////////////////////


#include "Wrappers.hh"
#include "FGeometryInit.hh"
#include "G4VPhysicalVolume.hh"
#include "FluggNavigator.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalVolumeStore.hh"
#include "globals.hh"

void lkfxwr(G4double& pSx, G4double& pSy, G4double& pSz,
            G4double* pV, const G4int& oldReg, const G4int& oldLttc,
	    G4int& newReg, G4int& flagErr, G4int& newLttc)
{
  //flag
#ifdef G4GEOMETRY_DEBUG
  G4cout << "======= LKFXWR =======" << G4endl;
#endif
  
  //FGeometryInit, navigator, volumeStore  pointers
  static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
  FluggNavigator * ptrNavig = ptrGeoInit->getNavigatorForTracking();
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  
  //coordinates in mm.
  G4ThreeVector pSource(pSx,pSy,pSz);
  pSource *= 10.0; //in millimeters!
  
  //locate point and update histories
  G4TouchableHistory * ptrTouchableHistory = 
    ptrGeoInit->GetTouchableHistory();
  ptrNavig->LocateGlobalPointAndUpdateTouchable(pSource,0,
						ptrTouchableHistory,true);
  //updating tmp but not old histories, they are useful in 
  //case of RGRPWR call, or when fluka, after a LOOKZ call, 
  //descards step for multiple scattering and returns to old history
  //NO, after lattice-fix we don't need old history anymore!
  
  ptrGeoInit->UpdateHistories(ptrTouchableHistory->GetHistory(),0); 
  G4VPhysicalVolume * located = ptrTouchableHistory->GetVolume();
  
  //if volume not found, out of mother volume: returns "number of volumes"+1
  if(!located) {
#ifdef G4GEOMETRY_DEBUG
    G4cout << "Out of mother volume!";
#endif
    G4int numVol = G4int(pVolStore->size());
    newReg = numVol + 1;
  }
  else { 
#ifdef G4GEOMETRY_DEBUG
    G4cout << "* ISVHWR call to store current NavHistWithCount in jrLtGeant"
	   << G4endl;
#endif 
    
    //save history in jrLtGeant and increment counter  
    G4int * jrLtGeant = ptrGeoInit->GetJrLtGeantArray();
    G4int LttcFlagGeant = ptrGeoInit->GetLttcFlagGeant();
    LttcFlagGeant += 1;
    jrLtGeant[LttcFlagGeant] = isvhwr(0,0);
    
#ifdef G4GEOMETRY_DEBUG
    G4cout << "* CONHWR call to increment counter" << G4endl;
#endif 
    G4int incrCount=1;
    conhwr(jrLtGeant[LttcFlagGeant],&incrCount);
    
    //update LttcFlagGeant
    ptrGeoInit->SetLttcFlagGeant(LttcFlagGeant);
    
    //return region number and dummy variables
    G4int volIndex = ~0;
    for (unsigned int i=0; i<pVolStore->size(); i++)
      if ((*pVolStore)[i] == located) 
	volIndex = i;
    // G4int volIndex=G4int(pVolStore->index(located));
    if (volIndex==(~0)) {
      G4cerr << "FLUGG: Problem in file WrapLookFX trying to find volume after step" << G4endl;
      exit(-999);
    }
    //G4int volIndex=G4int(pVolStore->index(located));
    newReg=volIndex+1;   
    newLttc=jrLtGeant[LttcFlagGeant];
    flagErr=newReg;
#ifdef G4GEOMETRY_DEBUG
    G4cout << "LKFXWR Located Physical volume = ";
    G4cout << located->GetName() << G4endl;
#endif
  }
}




