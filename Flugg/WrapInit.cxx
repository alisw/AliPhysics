
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// WrapInit.hh - Sara Vanini
//
// Wrapper for geometry initialisation.
//
// modified 12-IV-2000
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
//////////////////////////////////////////////////////////////////


#include "Wrappers.hh"
#include "FGeometryInit.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "globals.hh"

void jomiwr(const G4int & nge, const G4int& lin, const G4int& lou,
            G4int& flukaReg)
{
//flag
#ifdef G4GEOMETRY_DEBUG
  G4cout << "================== JOMIWR =================" << G4endl;
#endif 
  
  
  //Geoinit Pointer
#ifdef G4GEOMETRY_DEBUG
  G4cout << "\t *==> JOMIWR: Getting FGeometry..." << G4endl;
#endif 
  FGeometryInit * ptrGeoInit=FGeometryInit::GetInstance();
  
  //initialize geometry:construct detector and set world volume
#ifdef G4GEOMETRY_DEBUG
  G4cout << "\t *==> JOMIWR: Setting the detector..." << G4endl;
#endif 
//  ptrGeoInit->setDetector();
#ifdef G4GEOMETRY_DEBUG
  G4cout << "\t *==> JOMIWR: Setting mother volume..." << G4endl;
#endif 
//  ptrGeoInit->setMotherVolume(); 
  
  //close geometry for optimization
#ifdef G4GEOMETRY_DEBUG
  G4cout << "\t *==> JOMIWR: Closing geometry..." << G4endl;
#endif 
//  ptrGeoInit->closeGeometry();
  
  //initialize wrappers utility histories at the beginning of run and set flag
#ifdef G4GEOMETRY_DEBUG
  G4cout << "\t *==> JOMIWR: InitHistories..." << G4endl;
#endif 
//  ptrGeoInit->InitHistories();
  
  //initialize lattice array
#ifdef G4GEOMETRY_DEBUG
  G4cout << "\t *==> JOMIWR: Init lattice array..." << G4endl;
#endif 
//  ptrGeoInit->InitJrLtGeantArray();
  
  //initialize debug-array
#ifdef G4GEOMETRY_DEBUG
  G4cout << "\t *==> JOMIWR: Init debug array..." << G4endl;
#endif 
//  ptrGeoInit->InitHistArray();
  
  //create Fluka material cards in flukaMat.inp file
#ifdef G4GEOMETRY_DEBUG
  G4cout << "\t *==> JOMIWR: Init fluka materials..." << G4endl;
#endif 
//  ptrGeoInit->createFlukaMatFile();
  
  //returns number of volumes + 1
#ifdef G4GEOMETRY_DEBUG
  G4cout << "\t *==> JOMIWR: Returning..." << G4endl;
#endif 
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  G4int numVol = G4int(pVolStore->size());
  flukaReg = numVol + 1;
	
#ifdef G4GEOMETRY_DEBUG
  G4cout << "Number of volumes + 1: " << flukaReg << G4endl;
  G4cout << "================== Out of JOMIWR =================" << G4endl;
#endif
}



