// $Id$
// Flugg tag $Name$

#ifndef MyDetectorConstruction_h
#define MyDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();
     G4VPhysicalVolume* Construct();

  private:
     G4double expHall_rad;
     G4double expHall_z;

     G4double tar_rad;
     G4double tar_z;

     G4double litCil_rad;
     G4double litCil_z;

     G4double bigCil_rad;
     G4double bigCil_z;
};

#endif

