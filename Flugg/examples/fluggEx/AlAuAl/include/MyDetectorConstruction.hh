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
     G4double expHall_x;
     G4double expHall_y;
     G4double expHall_z;

     G4double myBox_x;
     G4double myBox_y;
     G4double myBox_zA;
     G4double myBox_zB;
     G4double myBox_zC;
     G4double myBox_zD;
     G4double myBox_zE;
};

#endif

