// $Id$
// Flugg tag $Name$

//
// MyDetectorConstruction.hh,  2000/02/11 for flugg
// Sara Vanini
//
// 


#ifndef MyDetectorConstruction_h
#define MyDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4Tubs;
class G4Sphere;
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class G4IntersectionSolid;


class MyDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    MyDetectorConstruction();
   ~MyDetectorConstruction();

  public:
     
     G4VPhysicalVolume* Construct();
  //     SetMagField(G4double fieldValue);

  private:
     
     G4Material*        TubMaterial;
     G4double           TubRad;
     
     G4Material*        SphereMaterial;
     G4double           SphereRad;
     
     G4int              NbOfLayers,NbOfTubs,NbOfSpheres;
     G4double           LayerThickness;
          
     G4double           DetSizeX,DetSizeY,DetSizeZ;
     
     G4Material*        defaultMaterial, *LayerMaterial;
     G4double           WorldSize;
            
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

     G4Box*             solidDet;    //pointer to the solid Det 
     G4LogicalVolume*   logicDet;    //pointer to the logical Det
     G4VPhysicalVolume* physiDet;    //pointer to the physical Det
     
     G4Box*             solidLayer;    //pointer to the solid Layer 
     G4LogicalVolume*   logicLayer;    //pointer to the logical Layer
     G4VPhysicalVolume* physiLayer1,* physiLayer2,* physiLayer3,
       * physiLayer4,* physiLayer5,* physiLayer6;    
       //pointer to the physical Layer
 
     G4Box*             solidTubLayer;    //pointer to the solid Layer 
     G4LogicalVolume*   logicTubLayer;    //pointer to the logical Layer
     G4VPhysicalVolume* physiTubLayer;    //pointer to the physical Layer
          
     G4Sphere*          solidSphere; //pointer to the solid Sphere
     G4LogicalVolume*   logicSphere; //pointer to the logical Sphere
     G4VPhysicalVolume* physiSphere; //pointer to the physical Sphere
     
     G4Tubs*            solidTub;      //pointer to the solid Tub
     G4LogicalVolume*   logicTub;      //pointer to the logical Tub
     G4VPhysicalVolume* physiTub;      //pointer to the physical Tub

     G4Box*             solidBoxIntersTub;
     G4VSolid*          solidTubSeg;      //pointer to the solid Tub Segment
     G4LogicalVolume*   logicTubSeg;      //pointer to the logical Tub Segment
     G4VPhysicalVolume* physiTubSeg;      //pointer to the physical Tub Segment

     G4Box*             solidBoxIntersSph;
     G4VSolid*          solidSphSeg1,*solidSphSeg2,*solidSphSeg3,*solidSphSeg4;
                                       //pointers to the solid Sphere Segment
     G4LogicalVolume*   logicSphSeg1,*logicSphSeg2,*logicSphSeg3,*logicSphSeg4;
                                       //pointer to the logical Sphere Segment
     G4VPhysicalVolume* physiSphSeg1,*physiSphSeg2,*physiSphSeg3,*physiSphSeg4;
                                       //pointer to the physical Sphere Segment



      
     G4UniformMagField* magField;      //pointer to the magnetic field
     
      
  private:
    
     void DefineMaterials();
     void ComputeDetectorParameters();
     G4VPhysicalVolume* ConstructDetector();     
};



inline void MyDetectorConstruction::ComputeDetectorParameters()
{
  // Compute derived parameters of the calorimeter
     LayerThickness = 2*SphereRad + 2*TubRad;
     DetSizeY       = LayerThickness * NbOfLayers;
}

#endif

