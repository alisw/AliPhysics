
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// WrapReg.hh - Sara Vanini 
//
// Wrapper for scoring hits: previous step end-point is taken from 
// history (and compared with fluka region index, flukaReg),
// then the wrapper returns all the information regarding the 
// volume tree, i.e. returns indMother[] array with all the 
// mother volumes index and repMother[] array with all the 
// mother volumes repetition number.   
//
// modified: 16/III/99
// modified: 14/IV/00 ptrLttc included 
// modified: 24.10.00: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
// modified: 17/06/02 by I. Gonzalez. STL migration.
//
///////////////////////////////////////////////////////////////////

#include "Wrappers.hh"
#include "FGeometryInit.hh"
#include "NavHistWithCount.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalVolumeStore.hh"
#include "globals.hh"


void rgrpwr(const G4int& flukaReg, const G4int& ptrLttc, G4int& g4Reg,
            G4int* indMother, G4int* repMother, G4int& depthFluka)
{
  //flag
#ifdef G4GEOMETRY_DEBUG
  G4cout << "============= RGRPWR ==============" << G4endl;    
  G4cout << "ptrLttc=" << ptrLttc << G4endl;
#endif 
  
  //Geoinit, Navigator, VolStore pointers
  static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  
  //get jrLtGeant array and flag
  G4int * jrLtGeant = ptrGeoInit->GetJrLtGeantArray();
  G4int LttcFlagGeant = ptrGeoInit->GetLttcFlagGeant();
  
  G4bool foundHistory = false;
  G4int i = LttcFlagGeant;
  
  while(!foundHistory && i>=0) {
    if(jrLtGeant[i]==ptrLttc) foundHistory = true;
    i -= 1;
  }
  
  if(!foundHistory) {
    G4cout << "* ERROR! History not in jrLtGeant!" << G4endl;
    //only in debugging version....
    assert(foundHistory);
  }
  else {
    //get history pointer from ptrLttc
    NavHistWithCount* ptrNavHistCount = 
      reinterpret_cast<NavHistWithCount*>(ptrLttc);
    G4NavigationHistory* ptrNavHist = ptrNavHistCount->GetNavHistPtr();
    
    G4VPhysicalVolume* ptrVolHistory = 0;
    G4int volHistIndex = 0;
    G4int depth = ptrNavHist->GetDepth();
    for(G4int h=0; h<=depth; h++) {
      ptrVolHistory = ptrNavHist->GetVolume(h);
      //
      volHistIndex = ~0;
      for (unsigned int i=0; i<pVolStore->size(); i++)
	if ((*pVolStore)[i] == ptrVolHistory) 
	  volHistIndex = i; 
      if (volHistIndex==(~0)) {
	G4cerr << "FLUGG: Problem in routine WrapReg tryingto find volume after step" << G4endl;
	exit(-999);
      }
      //volHistIndex = G4int(pVolStore->index(ptrVolHistory));
      //
      indMother[h] = volHistIndex+1;
      if(ptrVolHistory->IsReplicated()) {
	//true if volume is replica or parameterized,
	//false for placements; repetition numbers 
	//are set: 1,2,3,etc.; 0 for placed volumes.
	repMother[h] = 1+G4int(ptrNavHist->GetReplicaNo(h));
      }
      else {
	repMother[h] = 0;
      }
      
#ifdef  G4GEOMETRY_DEBUG
      //  G4cout<<"Level="<<h<<"   : ";
      //  G4cout<<"Region="<<indMother[h]<<"   (repetition="<<
      //    repMother[h]<<")"<<G4endl;
#endif
    }
    
    //compute new region index
    G4int volIndex = ~0;
    for (unsigned int i=0; i<pVolStore->size(); i++)
      if ((*pVolStore)[i] == ptrVolHistory) 
	volIndex = i;
    //G4int volIndex=G4int(pVolStore->index(ptrVolHistory));
    if (volIndex==(~0)) {
      G4cerr << "FLUGG: Problem in routine WrapReg tryingto find volume after step" << G4endl;
      exit(-999);
    }
    
    g4Reg=volIndex+1;
    
    depthFluka = depth;
  }
  
}





