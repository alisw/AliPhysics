// $Id$
// Flugg tag $Name$

#include "MyDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
  expHall_rad = 1000.*cm;
  expHall_z = 1000.*cm;

  tar_rad = 3.*cm;
  tar_z = 10.*cm;

  litCil_rad = 1000.*cm;
  litCil_z = 10.*cm;

  bigCil_rad = 1000.*cm;
  bigCil_z = 250.*cm;
}

MyDetectorConstruction::~MyDetectorConstruction()
{;}

G4VPhysicalVolume* MyDetectorConstruction::Construct()
{
//=======================volumes
//-------------------- experimental hall
  G4Tubs * experimentalHall_tub 
    = new G4Tubs("expHall_S",0.*cm,expHall_rad,expHall_z,0.*deg,360.*deg);
  G4LogicalVolume * experimentalHall_log
    = new G4LogicalVolume(experimentalHall_tub,0,"expHall_L",0,0,0);
  G4VPhysicalVolume * experimentalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"expHall_P",
			experimentalHall_log,0,false,0);

  //------------------------------ big cylinder
  G4Tubs * bigCil_tub
    = new G4Tubs("bigCil_S",0.*cm,bigCil_rad,bigCil_z,0.*deg,360.*deg);
  G4LogicalVolume * bigCil_log
    = new G4LogicalVolume(bigCil_tub,0,"bigCil_L",0,0,0);
  G4VPhysicalVolume * bigCil_phys
    = new G4PVPlacement(0,G4ThreeVector(0.*cm,0.*cm,750.*cm),"bigCil_P",
			bigCil_log,experimentalHall_phys,false,0);

  //------------------------------ target
  G4Tubs * target_tub
    = new G4Tubs("tar_S",0.*cm,tar_rad,tar_z,0.*deg,360.*deg);
  G4LogicalVolume * target_log
    = new G4LogicalVolume(target_tub,0,"tar_L",0,0,0);
  G4VPhysicalVolume * target_phys
    = new G4PVPlacement(0,G4ThreeVector(0.*cm,0.*cm,-10.*cm),"tar_P",
			target_log,experimentalHall_phys,false,0);

  //------------------------------ little cylinders
  G4Tubs * litCil_tub
    = new G4Tubs("litCil_S",0.*cm,litCil_rad,litCil_z,0.*deg,360.*deg);
  G4LogicalVolume * litCil_log
    = new G4LogicalVolume(litCil_tub,0,"litCil_L",0,0,0);

  /*  G4VPhysicalVolume * litCil_phy1
    = new G4PVPlacement(0,G4ThreeVector(0.*cm,0.*cm,10.*cm),
			"litCil_P1",litCil_log,experimentalHall_phys,false,0);
  G4VPhysicalVolume * litCil_phy2
    = new G4PVPlacement(0,G4ThreeVector(0.*cm,0.*cm,30.*cm),
			"litCil_P2",litCil_log,experimentalHall_phys,false,0);
  G4VPhysicalVolume * litCil_phy3
    = new G4PVPlacement(0,G4ThreeVector(0.*cm,0.*cm,50.*cm),
			"litCil_P3",litCil_log,experimentalHall_phys,false,0);
  G4VPhysicalVolume * litCil_phy4
    = new G4PVPlacement(0,G4ThreeVector(0.*cm,0.*cm,70.*cm),
			"litCil_P4",litCil_log,experimentalHall_phys,false,0);
			*/

       for(G4int i=0;i<25;i++)
  {
    cout<<"Cilindretto num. "<<i<<endl;
    new G4PVPlacement(0,G4ThreeVector(0.*cm,0.*cm,(1+2*i)*10.*cm),"litCil_P",
			litCil_log,experimentalHall_phys,false,i);
  }
  
  //------------------------------------------------------------------
  return experimentalHall_phys;
}
