// $Id$
// Flugg tag $Name$

#include "MyDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
  expHall_x = 10.*cm;
  expHall_y = 10.*cm;
  expHall_z = 10.*cm;

  myBox_x = 10.*cm;
  myBox_y = 10.*cm;
  myBox_zA = 5.*cm;
  myBox_zB = 0.00841665*cm;
  myBox_zC = 0.0010873*cm;
  myBox_zD = 0.07640375*cm;
  myBox_zE = 4.9140923*cm;

}

MyDetectorConstruction::~MyDetectorConstruction()
{;}

G4VPhysicalVolume* MyDetectorConstruction::Construct()
{
//===================================volumes
  G4cout << "MyDetectorConstruction::Construct start" << G4endl;


  //------------------------------ experimental hall
  G4Box * experimantalHall_box
    = new G4Box("expHall_b",expHall_x,expHall_y,expHall_z);
  G4LogicalVolume * experimantalHall_log
    = new G4LogicalVolume(experimantalHall_box,0,"expHall_L",0,0,0);
  G4VPhysicalVolume * experimantalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"expHall_P",
                        experimantalHall_log,0,false,0);

  //------------------------------ test box2
  G4Box * test_box2
    = new G4Box("myBox_b2",myBox_x,myBox_y,myBox_zA);
  G4LogicalVolume * test_log2
    = new G4LogicalVolume(test_box2,0,"myBox_L2",0,0,0);
  G4VPhysicalVolume * test_phys2
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,-5.*cm),"myBox_P2",
                        test_log2,experimantalHall_phys,false,0);

  //------------------------------ test box3
  G4Box * test_box3
    = new G4Box("myBox_b3",myBox_x,myBox_y,myBox_zB);
  G4LogicalVolume * test_log3
    = new G4LogicalVolume(test_box3,0,"myBox_L3",0,0,0);
  G4VPhysicalVolume * test_phys3
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.00841665*cm),"myBox_P3",
                        test_log3,experimantalHall_phys,false,0);

  //------------------------------ test box4
  G4Box * test_box4
    = new G4Box("myBox_b4",myBox_x,myBox_y,myBox_zC);
  G4LogicalVolume * test_log4
    = new G4LogicalVolume(test_box4,0,"myBox_L4",0,0,0);
  G4VPhysicalVolume * test_phys4
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.0179206*cm),"myBox_P4",
                        test_log4,experimantalHall_phys,false,0);

  //------------------------------ test box5
  G4Box * test_box5
    = new G4Box("myBox_b5",myBox_x,myBox_y,myBox_zD);
  G4LogicalVolume * test_log5
    = new G4LogicalVolume(test_box5,0,"myBox_L5",0,0,0);
  G4VPhysicalVolume * test_phys5
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.09541165*cm),"myBox_P5",
                        test_log5,experimantalHall_phys,false,0);

  //------------------------------ test box6
  G4Box * test_box6
    = new G4Box("myBox_b6",myBox_x,myBox_y,myBox_zE);
  G4LogicalVolume * test_log6
    = new G4LogicalVolume(test_box6,0,"myBox_L6",0,0,0);
  G4VPhysicalVolume * test_phys6
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,5.0859077*cm),"myBox_P6",
                        test_log6,experimantalHall_phys,false,0);  			

  //------------------------------------------------------------------
     G4int numLVVol = G4int(G4LogicalVolumeStore::GetInstance()->size());
     
  if (test_log6->GetMaterial())
    G4cout << "material is defined " << G4endl;
  else    
    G4cout << "material is NOT defined " << G4endl;

  return experimantalHall_phys;
}
