// $Id$
// Flugg tag $Name$

#include "FGeometryInit.hh"
#include "MyDetectorConstruction.hh"

#define flukam flukam_

extern "C" void flukam(const G4int & GeoFlag);

int main() {

  FGeometryInit* theFGeometryInit = FGeometryInit::GetInstance();
  
  // Test cases:
  // MyDetectorConstruction(Shape shape, G4int number, G4bool rotate);
  //
  // Shapes: 
  // kBox, kTubs,kCons, kTorus, kTrd, kTrap, kPara, kSphere,
  // kPolyhedra, kPolycone, kEllipticalTube, kHype  
  //
  // Number: in the interval < 1, 5 > 
  // =1: one volume without daughters 
  // >1: a volume with embedded daughters
  
  MyDetectorConstruction* detector 
    = new MyDetectorConstruction(kPara, 5, true);

  theFGeometryInit->setDetConstruction(detector);

//flag for geometry:
// 1 for GEANT4
// 0 for FLUKA
// 2 for Rubia
    const G4int flag = 1;

//call fortran
    flukam(flag);

//end
  return 0;
}




