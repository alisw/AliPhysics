// $Id$
// Flugg tag $Name$

#include "MyDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Para.hh"
#include "G4Sphere.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4EllipticalTube.hh"
#include "G4Hype.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"

#include "globals.hh"

//_____________________________________________________________________________
MyDetectorConstruction::MyDetectorConstruction(Shape shape, 
                                               G4int number,
					       G4bool rotate)
  : fAir(0),
    fWorld(0),

    fShape(shape),
    fNumber(number),
    fRotate(rotate),
    
    fWorldX(1.*m),
    fWorldY(1.*m),
    fWorldZ(1.*m),

    fX(5.*cm),
    fY(10.*cm),
    fZ(15.*cm),

    fR (5.*cm),
    fDR(5.*cm),
    fSphi (0.*deg),
    fDphi (315.*deg),
    fFullPhi(360.*deg),
  
    fAlpha(45.*deg),
    fTheta(30.*deg),
    fPhi  (15.*deg),

    fNofZPlanes(4),
    fNofSides(8)
{
  fRotation = new G4RotationMatrix();
  fRotation->rotateY(5.*deg);
}

//_____________________________________________________________________________
MyDetectorConstruction::~MyDetectorConstruction()
{;}

//
// private methods
//

//_____________________________________________________________________________
void MyDetectorConstruction::CreateAir()
{
  G4double a, iz, z, density;
  G4String name, symbol;
  G4double temperature, pressure;
  G4int nel;

  //Air
  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);
  density = 1.29*mg/cm3;
  fAir = new G4Material(name="Air", density, nel=2);
  fAir->AddElement(elN, .7);
  fAir->AddElement(elO, .3);
}


//_____________________________________________________________________________
void MyDetectorConstruction::CreateWorld()
{
  
  G4Box* box 
    = new G4Box("world_s", fWorldX, fWorldY, fWorldZ);
    
  G4LogicalVolume* log 
    = new G4LogicalVolume(box, fAir, "world_l", 0, 0, 0);
    
  fWorld
    = new G4PVPlacement(0, G4ThreeVector(), "world_p", log, 0, false, 0);

}

//_____________________________________________________________________________
void MyDetectorConstruction::CreateBox()
{
  
  // big box

  G4Box* solid_env 
    = new G4Box("box_env_s", fX, fY, fZ*fNumber);
    
  G4LogicalVolume* log_env 
    = new G4LogicalVolume(solid_env, fAir, "box_env_l", 0, 0, 0);
    
  // rotation        
  G4RotationMatrix* rotation = 0;
  if (fRotate) rotation = fRotation;

  new G4PVPlacement(rotation, G4ThreeVector(0.,0.,0.),
                    "box_env_p", log_env, fWorld, false, 0);

  // small daughter boxes

  if (fNumber>1) {
     G4Box* box 
       = new G4Box("box_s", fX, fY, fZ);
    
     G4LogicalVolume* log 
       = new G4LogicalVolume(box, fAir, "box_l", 0, 0, 0);
    
    for (G4int i=0; i<fNumber; i++) {
  
      // position
      G4double zpos = -fNumber*fZ + fZ + i*2*fZ;
      G4ThreeVector position(0., 0., zpos);
    
      new G4PVPlacement(0, position, log, "box_p", log_env, false, 0);
    }  
  }
}

//_____________________________________________________________________________
void MyDetectorConstruction::CreateTubs()
{
  
  // big tub segment

  G4Tubs* solid_env 
    = new G4Tubs("tubs_env_s", fR, fR+fDR*fNumber, fZ, fSphi, fDphi);
    
  G4LogicalVolume* log_env 
    = new G4LogicalVolume(solid_env, fAir, "tubs_env_l", 0, 0, 0);
    
  // rotation        
  G4RotationMatrix* rotation = 0;
  if (fRotate) rotation = fRotation;

  new G4PVPlacement(rotation, G4ThreeVector(0.,0.,0.),
                    "tubs_env_p", log_env, fWorld, false, 0);

  // small daughter tub segments

  if (fNumber>1) {
    for (G4int i=0; i<fNumber; i++) {
      G4Tubs* solid 
        = new G4Tubs("tubs_s", fR+i*fDR, fR+(i+1)*fDR, fZ, fSphi, fDphi);
    
      G4LogicalVolume* log 
        = new G4LogicalVolume(solid, fAir, "tubs_l", 0, 0, 0);
    
      // position
      G4ThreeVector position(0., 0., 0.);
    
      new G4PVPlacement(0, position, log, "tubs_p", log_env, false, 0);
    }  
  }
}

//_____________________________________________________________________________
void MyDetectorConstruction::CreateTorus()
{
  
  // big tourus segment

  G4Torus* solid_env 
    = new G4Torus("torus_env_s", 
                  fR, fR+fDR, fR+5*fDR, fSphi, fDphi);
    
  G4LogicalVolume* log_env 
    = new G4LogicalVolume(solid_env, fAir, "torus_env_l", 0, 0, 0);
    
  // rotation        
  G4RotationMatrix* rotation = 0;
  if (fRotate) rotation = fRotation;

  new G4PVPlacement(rotation, G4ThreeVector(0.,0.,0.),
                    "torus_env_p", log_env, fWorld, false, 0);

  // small daughter torus segments
  if (fNumber>1) {
    for (G4int i=0; i<fNumber; i++) {
      G4Torus* solid 
        = new G4Torus("torus_s",
	              fR, fR+fDR, fR+5*fDR, 
		      fSphi+i*(fDphi/fNumber), fDphi/fNumber);
    
      G4LogicalVolume* log 
        = new G4LogicalVolume(solid, fAir, "torus_l", 0, 0, 0);
    
      // position
      G4ThreeVector position(0., 0., 0.);
    
      new G4PVPlacement(0, position, log, "torus_p", log_env, false, 0);
    }  
  }
}

//_____________________________________________________________________________
void MyDetectorConstruction::CreateCons()
{
  
  // big cons segment

  G4Cons* solid_env 
    = new G4Cons("cons_env_s", 
                 fR, fR+fDR*fNumber, fR+fDR, fR+fDR+fDR*fNumber, 
		 fZ, fSphi, fDphi);
    
  G4LogicalVolume* log_env 
    = new G4LogicalVolume(solid_env, fAir, "cons_env_l", 0, 0, 0);
    
  // rotation        
  G4RotationMatrix* rotation = 0;
  if (fRotate) rotation = fRotation;

  new G4PVPlacement(rotation, G4ThreeVector(0.,0.,0.),
                    "cons_env_p", log_env, fWorld, false, 0);

  // small daughter cons segments

  if (fNumber>1) {
    for (G4int i=0; i<fNumber; i++) {
      G4Cons* solid 
        = new G4Cons("cons_s", 
	             fR+i*fDR, fR+(i+1)*fDR, fR+fDR+i*fDR, fR+fDR+(i+1)*fDR,
		     fZ, fSphi, fDphi);
    
      G4LogicalVolume* log 
        = new G4LogicalVolume(solid, fAir, "cons_l", 0, 0, 0);
    
      // position
      G4ThreeVector position(0., 0., 0.);
    
      new G4PVPlacement(0, position, log, "cons_p", log_env, false, 0);
    }  
  }
}

//_____________________________________________________________________________
void MyDetectorConstruction::CreateTrd()
{
  
  // big trd segment

  G4Trd* solid_env 
    = new G4Trd("trd_env_s", 2*fX, 4*fX, 2*fY, 4*fY, fZ*fNumber);
    
  G4LogicalVolume* log_env 
    = new G4LogicalVolume(solid_env, fAir, "trd_env_l", 0, 0, 0);
    
  // rotation        
  G4RotationMatrix* rotation = 0;
  if (fRotate) rotation = fRotation;

  new G4PVPlacement(rotation, G4ThreeVector(0.,0.,0.),
                    "trd_env_p", log_env, fWorld, false, 0);

  // small daughter trd segments

  G4Trd* solid 
    = new G4Trd("trd_s", fX, 2*fX, fY, 2*fY, fZ);

  G4LogicalVolume* log 
    = new G4LogicalVolume(solid, fAir, "trd_l", 0, 0, 0);

  if (fNumber>1) {
    for (G4int i=0; i<fNumber; i++) {
      // position
      G4double zpos = -fNumber*fZ + fZ + i*2*fZ;
      G4ThreeVector position(0., 0., zpos);
    
      new G4PVPlacement(0, position, log, "trd_p", log_env, false, 0);
    }  
  }
}

//_____________________________________________________________________________
void MyDetectorConstruction::CreateTrap()
{
  // big trap

  G4double theta = fTheta;
  G4double phi = fPhi;
  if (fNumber>1) {
    theta = 0.;
    phi = 0.;
  }

  G4Trap* solid_env 
    = new G4Trap("trap_env_s", 
                 fZ*fNumber, theta, phi,
		 fY, fX, 2*fX, 5.0*deg, fY, fX, 2*fX, 5.0*deg);
    
  G4LogicalVolume* log_env 
    = new G4LogicalVolume(solid_env, fAir, "trap_env_l", 0, 0, 0);
    
  // rotation        
  G4RotationMatrix* rotation = 0;
  if (fRotate) rotation = fRotation;

  new G4PVPlacement(rotation, G4ThreeVector(0.,0.,0.),
                    "trap_env_p", log_env, fWorld, false, 0);

  // small daughter trap segments
  
  G4Trap* solid 
    = new G4Trap("trap_env_s", 
                 fZ,  theta, phi,
		 fY, fX, 2*fX, 5.0*deg, fY, fX, 2*fX, 5.0*deg);

  G4LogicalVolume* log 
    = new G4LogicalVolume(solid, fAir, "trap_l", 0, 0, 0);

  if (fNumber>1) {
    for (G4int i=0; i<fNumber; i++) {
      // position
      G4double zpos = -fNumber*fZ + fZ + i*2*fZ;
      G4ThreeVector position(0., 0., zpos);
    
      new G4PVPlacement(0, position, log, "trap_p", log_env, false, 0);
    }  
  }
}

//_____________________________________________________________________________
void MyDetectorConstruction::CreatePara()
{
  
  G4double theta = fTheta;
  G4double phi = fPhi;
  if (fNumber>1) {
    theta = 0.;
    phi = 0.;
  }

  // big para segment

  G4Para* solid_env 
    = new G4Para("para_env_s", fX, fY, fZ*fNumber, fAlpha, theta, phi);
    
  G4LogicalVolume* log_env 
    = new G4LogicalVolume(solid_env, fAir, "para_env_l", 0, 0, 0);
    
  // rotation        
  G4RotationMatrix* rotation = 0;
  if (fRotate) rotation = fRotation;

  new G4PVPlacement(rotation, G4ThreeVector(0.,0.,0.),
                    "para_env_p", log_env, fWorld, false, 0);

  // small daughter para segments

  G4Para* solid 
    = new G4Para("para_s", fX, fY, fZ, fAlpha, theta, phi);

  G4LogicalVolume* log 
    = new G4LogicalVolume(solid, fAir, "para_l", 0, 0, 0);

  if (fNumber>1) {
    for (G4int i=0; i<fNumber; i++) {
      // position
      G4double zpos = -fNumber*fZ + fZ + i*2*fZ;
      G4ThreeVector position(0., 0., zpos);
    
      new G4PVPlacement(0, position, log, "para_p", log_env, false, 0);
    }  
  }  
}

//_____________________________________________________________________________
void MyDetectorConstruction::CreateSphere()
{
  
  // big sphere segment

  G4Sphere* solid_env 
    = new G4Sphere("sphere_env_s", 
                 fR, fR+fDR*fNumber, fSphi, fDphi, fSphi, fDphi/2.); 
    
  G4LogicalVolume* log_env 
    = new G4LogicalVolume(solid_env, fAir, "sphere_env_l", 0, 0, 0);
    
  // rotation        
  G4RotationMatrix* rotation = 0;
  if (fRotate) rotation = fRotation;

  new G4PVPlacement(rotation, G4ThreeVector(0.,0.,0.),
                    "sphere_env_p", log_env, fWorld, false, 0);

  // small daughter sphere segments
  if (fNumber>1) {
    for (G4int i=0; i<fNumber; i++) {
      G4Sphere* solid 
        = new G4Sphere("sphere_s", 
	             fR+i*fDR, fR+(i+1)*fDR, fSphi, fDphi, fSphi, fDphi/2.);
    
      G4LogicalVolume* log 
        = new G4LogicalVolume(solid, fAir, "sphere_l", 0, 0, 0);
    
      // position
      G4ThreeVector position(0., 0., 0.);
    
      new G4PVPlacement(0, position, log, "sphere_p", log_env, false, 0);
    }  
  }  
}

//_____________________________________________________________________________
void MyDetectorConstruction::CreatePolyhedra()
{
  // big polyhedra

  G4double* z_env    = new G4double[fNofZPlanes];
  G4double* rmin_env = new G4double[fNofZPlanes];
  G4double* rmax_env = new G4double[fNofZPlanes];
  
  for (G4int i=0; i< fNofZPlanes; i++) {
    z_env[i]    = -fNofZPlanes*fZ + fZ + i*2*fZ;
    rmin_env[i] = fR + (i+1) * fDR + (2*(i%2)-1) * (i+1) * fDR/2.;
    rmax_env[i] = fR + (i+2) * fDR + (2*(i%2)-1) * (i+1) * fDR/2.;
  }  

  G4Polyhedra* solid_env 
    = new G4Polyhedra("polyhedra_env_s", fSphi, fFullPhi,
                      fNofSides, fNofZPlanes, z_env, rmin_env, rmax_env);
    
  G4LogicalVolume* log_env 
    = new G4LogicalVolume(solid_env, fAir, "polyhedra_env_l", 0, 0, 0);

  // rotation        
  G4RotationMatrix* rotation = 0;
  if (fRotate) rotation = fRotation;

  new G4PVPlacement(rotation, G4ThreeVector(0.,0.,0.),
                    "polyhedra_env_p", log_env, fWorld, false, 0);

  // small daughter polyhedra segments

  if (fNumber>1) {
    for (G4int j=0; j<fNumber; j++) {
    
      G4double* z    = new G4double[fNofZPlanes];
      G4double* rmin = new G4double[fNofZPlanes];
      G4double* rmax = new G4double[fNofZPlanes];
  
      for (G4int k=0; k< fNofZPlanes; k++) {
	G4double dr = (rmax_env[k] - rmin_env[k])/fNumber;
        z[k]    = z_env[k];
        rmin[k] = rmin_env[k] + j*dr;
        rmax[k] = rmin_env[k] + (j+1)*dr;
      }  
    
      G4Polyhedra* solid 
        = new G4Polyhedra("polyhedra_s", fSphi, fFullPhi,
                          fNofSides, fNofZPlanes, z, rmin, rmax);
    
      G4LogicalVolume* log 
        = new G4LogicalVolume(solid, fAir, "polyhedra_l", 0, 0, 0);
    
      // position
      G4ThreeVector position(0., 0., 0.);
      new G4PVPlacement(0, position, log, "polyhedra_p", log_env, false, 0);
      
      delete [] z;
      delete [] rmin;
      delete [] rmax;
    }  
  }

  delete [] z_env;
  delete [] rmin_env;
  delete [] rmax_env;
}

//_____________________________________________________________________________
void MyDetectorConstruction::CreatePolycone()
{
  // big polycone

  G4double* z_env    = new G4double[fNofZPlanes];
  G4double* rmin_env = new G4double[fNofZPlanes];
  G4double* rmax_env = new G4double[fNofZPlanes];
  
  for (G4int i=0; i< fNofZPlanes; i++) {
    z_env[i]    = -fNofZPlanes*fZ + fZ + i*2*fZ;
    rmin_env[i] = fR + (i+1) * fDR + (2*(i%2)-1) * (i+1) * fDR/2.;
    rmax_env[i] = fR + (i+2) * fDR + (2*(i%2)-1) * (i+1) * fDR/2.;
  }  

  G4Polycone* solid_env 
    = new G4Polycone("polycone_env_s", fSphi, fFullPhi,
                      fNofZPlanes, z_env, rmin_env, rmax_env);
    
  G4LogicalVolume* log_env 
    = new G4LogicalVolume(solid_env, fAir, "polycone_env_l", 0, 0, 0);

  // rotation        
  G4RotationMatrix* rotation = 0;
  if (fRotate) rotation = fRotation;

  new G4PVPlacement(rotation, G4ThreeVector(0.,0.,0.),
                    "polycone_env_p", log_env, fWorld, false, 0);

  // small daughter polycone segments

  if (fNumber>1) {
    for (G4int j=0; j<fNumber; j++) {
    
      G4double* z    = new G4double[fNofZPlanes];
      G4double* rmin = new G4double[fNofZPlanes];
      G4double* rmax = new G4double[fNofZPlanes];
  
      for (G4int k=0; k< fNofZPlanes; k++) {
	G4double dr = (rmax_env[k] - rmin_env[k])/fNumber;
        z[k]    = z_env[k];
        rmin[k] = rmin_env[k] + j*dr;
        rmax[k] = rmin_env[k] + (j+1)*dr;
      }  
    
      G4Polycone* solid 
        = new G4Polycone("polycone_s", fSphi, fFullPhi,
                          fNofZPlanes, z, rmin, rmax);
    
      G4LogicalVolume* log 
        = new G4LogicalVolume(solid, fAir, "polycone_l", 0, 0, 0);
    
      // position
      G4ThreeVector position(0., 0., 0.);
      new G4PVPlacement(0, position, log, "polycone_p", log_env, false, 0);
      
      delete [] z;
      delete [] rmin;
      delete [] rmax;
    }  
  }

  delete [] z_env;
  delete [] rmin_env;
  delete [] rmax_env;
}

//_____________________________________________________________________________
void MyDetectorConstruction::CreateEllipticalTube()
{
  
  // big elliptical tube
  // visualization not available

  G4EllipticalTube* solid_env 
    = new G4EllipticalTube("eltube_env_s", 
                           fR + 2.*fDR/5., fR + 4.*fDR/5., fZ);
    
  G4LogicalVolume* log_env 
    = new G4LogicalVolume(solid_env, fAir, "eltube_env_l", 0, 0, 0);
    
  // rotation        
  G4RotationMatrix* rotation = 0;
  if (fRotate) rotation = fRotation;

  new G4PVPlacement(rotation, G4ThreeVector(0.,0.,0.),
                    "eltube_env_p", log_env, fWorld, false, 0);

  // smaller daughter elliptical tube
  if (fNumber>1) {
    G4EllipticalTube* solid 
      = new G4EllipticalTube("eltube_s", fR, fR + fDR/5., fZ);
    
    G4LogicalVolume* log 
      = new G4LogicalVolume(solid, fAir, "eltube_l", 0, 0, 0);
    
    // position
    G4ThreeVector position(0., 0., 0.);
    
    new G4PVPlacement(0, position, log, "eltube_p", log_env, false, 0);
  }

  if (fNumber>2) {
    // not supported cases
    G4cout << "Elliptical tube:" << G4endl;
    G4cout << "Only one daugher is implemented in the test." << G4endl;
  }  
}

//_____________________________________________________________________________
void MyDetectorConstruction::CreateHype()
{
  
  // hyperbolic volume
  // visualization not available 

  G4Hype* solid_env 
    = new G4Hype("hype_env_s", 
                 fR, fR + fDR, fPhi, fTheta, fZ);
    
  G4LogicalVolume* log_env 
    = new G4LogicalVolume(solid_env, fAir, "hype_env_l", 0, 0, 0);
    
  // rotation        
  G4RotationMatrix* rotation = 0;
  if (fRotate) rotation = fRotation;

  new G4PVPlacement(rotation, G4ThreeVector(0.,0.,0.),
                    "hype_env_p", log_env, fWorld, false, 0);
  if (fNumber>1) {
    // not supported cases
    G4cout << "Hype:" << G4endl;
    G4cout << "No daughers are implemented in the test." << G4endl;
  } 
}

//
// public methods
//

//_____________________________________________________________________________
G4VPhysicalVolume* MyDetectorConstruction::Construct()
{
  CreateAir();

  CreateWorld();
  
  switch (fShape) {

    // CSG solids
    case kBox:
           CreateBox();
	   break;
    case kTubs:
           CreateTubs();
	   break;
    case kCons:
           CreateCons();
	   break;
    case kTorus:
           CreateTorus();
	   break;
    case kTrd:
           CreateTrd();
	   break;
    case kTrap:
           CreateTrap();
	   break;
    case kPara:
           CreatePara();
	   break;
    case kSphere:
           CreateSphere();
	   break;

    // specific solids
    case kPolyhedra:
           CreatePolyhedra();
	   break;
    case kPolycone:
           CreatePolycone();
	   break;
    case kEllipticalTube:
           CreateEllipticalTube();
	   break;
    case kHype:
           CreateHype();
	   break;
  }	   

  return fWorld;
}
