// $Id$
// Flugg tag $Name$

// MyDetectorConstruction.hh, 11/XII/1998  -  Sara Vanini
// GEANT4 tag $Name$
//
// 

#ifndef MyDetectorConstruction_h
#define MyDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
//class T36EMCalorimeterSD;
//class T36HADCalorimeterSD;

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    MyDetectorConstruction();
   ~MyDetectorConstruction();
    G4VPhysicalVolume* Construct();
    void SetMagField(G4double);
    void PrintCalorParameters();

//for N02 (1+staies for the last HAD absorber):
    G4int GetNbOfLayers()      {return (1+NbOfEMLayers+NbOfHADLayers);};


    G4int GetNbOfEMLayers()              {return NbOfEMLayers;};
    G4int GetNbOfHADLayers()              {return NbOfHADLayers;};
    G4double GetWorldSizeX()           {return WorldSizeX;};
    G4double GetWorldSizeY()          {return WorldSizeY;};
    G4double GetWorldSizeZ()          {return WorldSizeZ;};
    G4double GetCalorSizeY()          {return ModuleSizeY;};
    G4double GetCalorSizeZ()          {return CalorSizeZ;};
    const G4VPhysicalVolume* GetEMphysiAbs()   {return EMphysiAbsorber;};
    const G4VPhysicalVolume* GetEMphysiGap()        {return EMphysiGap;};
    const G4VPhysicalVolume* GetHADphysiAbs()   {return HADphysiAbsorber;};
    const G4VPhysicalVolume* GetlastHADAbsorber()   {return lastHADphysiAbsorber;};
    const G4VPhysicalVolume* GetHADphysiGap()        {return HADphysiGap;};
        
    const G4VPhysicalVolume* GetEMphysiAbMedShield() {return EMphysiAbMedShield;};

    const G4VPhysicalVolume* GetEMphysiSc() {return EMphysiSc;};

    const G4VPhysicalVolume* GetHADphysiAbMedShield() {return HADphysiAbMedShield;};
    const G4VPhysicalVolume* GetHADphysiAbExShield() {return HADphysiAbExShield;};
    const G4VPhysicalVolume* GetHADphysiSc() {return HADphysiSc;};

    const G4VPhysicalVolume* GetlastHADphysiAbs() {return lastHADphysiAbsorber;};
   
        
  private:
     
     G4double           EMModuleThickness,EMModuleSizeZ;
     G4double           HADModuleThickness,HADModuleSizeZ;
     G4Material*        AbsorberMaterial,*GapMaterial,*defaultMaterial;
     G4Material*        RodMaterial,*ScMaterial,*SpacerMaterial;
     G4Material*        WLSMaterial,*frontPlateMaterial,*MedAbsMaterial,*ExAbsMaterial; 
     G4Material*        vacuum,*WLSMat,*Al,*MedLead,*ExLead;
     G4double           AbsorberThickness,AbMedShieldThickness,AbExShieldThickness;
     
     G4double           GapThickness;
     
     G4int              NbOfEMLayers,NbOfHADLayers,NbOfModules;
     G4double           LayerThickness;
          
     G4double           ModuleSizeX,ModuleSizeY,ModuleSizeZ,CalorSizeZ;
     G4double		EMWLSSizeX,EMWLSSizeZ,RodRad;
     G4double		HADWLSSizeX,HADWLSSizeZ;
     G4double		ScThickness,ScSizeY,AlSizeX,SpaSizeX;
     G4double		SpacerThickness,SpacerSizeY;

     G4double           WorldSizeX,WorldSizeY,WorldSizeZ;
            
     G4Box*             solidWorld,*solidFrontVac;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld,*logicFrontVac;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld,*physiFrontVac;    //pointer to the physical World

     G4Box*             solidCalor,*solidModule,*solidEMWLS,*solidHADWLS,*solidSpace; 
     G4LogicalVolume*   logicCalor,*logicModule,*logicEMWLS,*logicHADWLS,*logicSpace;
     G4VPhysicalVolume* physiCalor,*physiModule,*physiSpace;
     
     G4Box*             EMsolidLayer,*HADsolidLayer,*solidAl,*EMsolidModule,*HADsolidModule; 
     G4LogicalVolume*   EMlogicLayer,*HADlogicLayer,*logicAl,*EMlogicModule,*HADlogicModule;
     G4VPhysicalVolume* EMphysiLayer,*HADphysiLayer,*physiAl,*EMphysiModule,*HADphysiModule;
         
     G4Box*             EMsolidAbsorber,*HADsolidAbsorber,*lastHADsolidAbsorber;
     G4LogicalVolume*   EMlogicAbsorber,*HADlogicAbsorber,*lastHADlogicAbsorber; 
     G4VPhysicalVolume* EMphysiAbsorber,*HADphysiAbsorber,*lastHADphysiAbsorber;
     
     G4Box*             EMsolidAbMedShield,*HADsolidAbMedShield,*EMsolidAbExShield,*HADsolidAbExShield;
     G4LogicalVolume*   EMlogicAbMedShield,*HADlogicAbMedShield,*EMlogicAbExShield,*HADlogicAbExShield; 
     G4VPhysicalVolume* EMphysiAbMedShield,*HADphysiAbMedShield,*lastHADphysiAbMedShield;
     G4VPhysicalVolume* EMphysiAbExShield,*HADphysiAbExShield,*lastHADphysiAbExShield;   

     G4Box*             EMsolidGap,*HADsolidGap,*EMsolidSc,*HADsolidSc;
     G4LogicalVolume*   EMlogicGap,*HADlogicGap,*EMlogicSc,*HADlogicSc;
     G4VPhysicalVolume* EMphysiGap,*HADphysiGap,*EMphysiSc,*HADphysiSc;

     G4Box* 		EMsolidSpacer,*HADsolidSpacer;
     G4Tubs*		EMsolidRod,*HADsolidRod;
     G4LogicalVolume*   EMlogicRod,*HADlogicRod,*EMlogicSpacer,*HADlogicSpacer;
     
     G4UniformMagField* magField;      //pointer to the magnetic field

     //pointer to the sensitive detectors - EM and HAD
  //     T36EMCalorimeterSD* EMcalorimeterSD;
  //     T36HADCalorimeterSD* HADcalorimeterSD;  

//for N02
//     ExN02CalorimeterSD* calorimeterSD;
 
  private:
   
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};



inline void MyDetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter

     LayerThickness = AbsorberThickness + GapThickness;
     EMModuleThickness = NbOfEMLayers*LayerThickness;
     HADModuleThickness = NbOfHADLayers*LayerThickness;

     EMWLSSizeX = EMModuleThickness + HADModuleThickness + AbsorberThickness;
     HADWLSSizeX = HADModuleThickness + AbsorberThickness;

     ModuleSizeZ = 2*EMWLSSizeZ + EMModuleSizeZ;
     CalorSizeZ = NbOfModules*ModuleSizeZ;
     ModuleSizeX = AlSizeX + (NbOfEMLayers+NbOfHADLayers)*LayerThickness
		 + AbsorberThickness + SpaSizeX; 
}

#endif

