// $Id$
// Flugg tag $Name$

#ifndef MyDetectorConstruction_h
#define MyDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4Material;

enum Shape {

  // CSG solids 
  kBox, kTubs,kCons, kTorus, kTrd, kTrap, kPara, kSphere,
       
  // specific solids
  kPolyhedra, kPolycone, kEllipticalTube, kHype 
};


class MyDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    MyDetectorConstruction(Shape shape, G4int number,G4bool rotate);
    virtual ~MyDetectorConstruction();

    // methods
    virtual G4VPhysicalVolume* Construct();

  private:
    MyDetectorConstruction() {}
  
    // methods
    void CreateAir();
    void CreateWorld();
    
    // CSG solids
    void CreateBox();
    void CreateTubs();
    void CreateCons();
    void CreateTorus();
    void CreateTrd();
    void CreateTrap();
    void CreatePara();
    void CreateSphere();
  
    // specific solids
    void CreatePolyhedra();
    void CreatePolycone();
    void CreateEllipticalTube();
    void CreateHype();

    // data members
    G4Material*        fAir;
    G4VPhysicalVolume* fWorld;
    G4RotationMatrix*  fRotation;

    Shape   fShape;
    G4int   fNumber;
    G4bool  fRotate;

    const G4double fWorldX;
    const G4double fWorldY;
    const G4double fWorldZ;

    const G4double fX;
    const G4double fY;
    const G4double fZ;

    const G4double fR;
    const G4double fDR;
    const G4double fSphi;
    const G4double fDphi;
    const G4double fFullPhi;

    const G4double fAlpha;
    const G4double fTheta;
    const G4double fPhi;

    const G4int fNofZPlanes;
    const G4int fNofSides;
};

#endif

