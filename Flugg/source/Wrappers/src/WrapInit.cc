// $Id$
// Flugg tag $Name$

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
	G4cout<<"================== JOMIWR ================="<<G4endl;
#endif 

       
//Geoinit Pointer
        FGeometryInit * ptrGeoInit=FGeometryInit::GetInstance();
	
//initialize geometry:construct detector and set world volume
	ptrGeoInit->setDetector();
	ptrGeoInit->setMotherVolume(); 

#ifdef G4GEOMETRY_DEBUG
	G4cout << "settings1 " << G4endl;
#endif

//close geometry for optimization
        ptrGeoInit->closeGeometry();

#ifdef G4GEOMETRY_DEBUG
	G4cout << "close geometry "<< flukaReg << G4endl;
#endif

//initialize wrappers utility histories at the beginning of run and set flag

	ptrGeoInit->InitHistories();

#ifdef G4GEOMETRY_DEBUG
	G4cout << " InitHistories" << G4endl;
#endif

//initialize lattice array
	ptrGeoInit->InitJrLtGeantArray();

#ifdef G4GEOMETRY_DEBUG
	G4cout << "InitJrLtGeantArray " << G4endl;
#endif

//initialize debug-array
	ptrGeoInit->InitHistArray();

#ifdef G4GEOMETRY_DEBUG
	G4cout << "InitHistArray " << G4endl;
#endif


//create Fluka material cards in flukaMat.inp file
       	ptrGeoInit->createFlukaMatFile();

#ifdef G4GEOMETRY_DEBUG
	G4cout<<"createFlukaMatFile "<<flukaReg<<G4endl;
#endif

//returns number of volumes + 1
        G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
	G4int numVol = G4int(pVolStore->size());
	flukaReg = numVol + 1;
	
#ifdef G4GEOMETRY_DEBUG
	G4cout<<"Number of volumes + 1: "<<flukaReg<<G4endl;
#endif
}



