// $Id$
// Flugg tag $Name$

//
// MyDetectorConstruction.hh, 3/III/1998  -  Sara Vanini
// T36 calorimeter with parametric volumes insted of replicans
//
// 

//#include "T36EMCalorimeterSD.hh"
//#include "T36HADCalorimeterSD.hh"
#include "MyDetectorConstruction.hh"
//#include "T36ModuleParam.hh"
#include "T36LayerParam.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
//#include "G4SDManager.hh"
//#include "G4RunManager.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


MyDetectorConstruction::MyDetectorConstruction()
:solidWorld(NULL),logicWorld(NULL),physiWorld(NULL),
 solidCalor(NULL),logicCalor(NULL),physiCalor(NULL),
 solidModule(NULL),logicModule(NULL),physiModule(NULL),
 EMsolidModule(NULL),EMlogicModule(NULL),EMphysiModule(NULL),
 HADsolidModule(NULL),HADlogicModule(NULL),HADphysiModule(NULL),
 solidFrontVac(NULL),logicFrontVac(NULL),physiFrontVac(NULL),

 solidEMWLS(NULL),logicEMWLS(NULL),solidHADWLS(NULL),logicHADWLS(NULL),
 solidAl(NULL),logicAl(NULL),physiAl(NULL),
 solidSpace(NULL),logicSpace(NULL),physiSpace(NULL),
 EMsolidLayer(NULL),EMlogicLayer(NULL),EMphysiLayer(NULL),
 HADsolidLayer(NULL),HADlogicLayer(NULL),HADphysiLayer(NULL),

 EMsolidAbsorber(NULL),EMlogicAbsorber(NULL),EMphysiAbsorber(NULL),
 HADsolidAbsorber(NULL),HADlogicAbsorber(NULL),HADphysiAbsorber(NULL),
 lastHADsolidAbsorber(NULL),lastHADlogicAbsorber(NULL),lastHADphysiAbsorber(NULL),
 lastHADphysiAbMedShield(NULL),lastHADphysiAbExShield(NULL),
 EMsolidAbMedShield(NULL),EMlogicAbMedShield(NULL),EMphysiAbMedShield(NULL),
 EMsolidAbExShield(NULL),EMlogicAbExShield(NULL),EMphysiAbExShield(NULL),
 HADsolidAbMedShield(NULL),HADlogicAbMedShield(NULL),HADphysiAbMedShield(NULL),
 HADsolidAbExShield(NULL),HADlogicAbExShield(NULL),HADphysiAbExShield(NULL),


 EMsolidGap(NULL),EMlogicGap(NULL),EMphysiGap(NULL),
 HADsolidGap(NULL),HADlogicGap(NULL),HADphysiGap(NULL),
 EMsolidSc(NULL),EMlogicSc(NULL),EMphysiSc(NULL),
 HADsolidSc(NULL),HADlogicSc(NULL),HADphysiSc(NULL),
 AbsorberMaterial(NULL),GapMaterial(NULL),ScMaterial(NULL),WLSMaterial(NULL),
 RodMaterial(NULL),defaultMaterial(NULL),SpacerMaterial(NULL),
 frontPlateMaterial(NULL),MedAbsMaterial(NULL),ExAbsMaterial(NULL),
 EMsolidRod(NULL),HADsolidRod(NULL),EMlogicRod(NULL),HADlogicRod(NULL),
 EMsolidSpacer(NULL),HADsolidSpacer(NULL),EMlogicSpacer(NULL),HADlogicSpacer(NULL)

{
  WorldSizeX = 4000*mm; 
  WorldSizeY = 4000*mm;
  WorldSizeZ = 4000*mm;    

  // default parameter values of the calorimeter
  // EMWLSSizeX       = 1103.5*mm; //copre tutta la torre
  EMWLSSizeZ       = 5*mm;    
  // HADWLSSizeX      = 887.5*mm;  //come input fluka: dal layer 17 per tutta la torre
  HADWLSSizeZ      = 3.5*mm;  
  ScSizeY          = 218*mm;
  ScThickness      = 2.5*mm;
  RodRad           = 1*mm;
  AlSizeX          = 20*mm;
  SpaSizeX         = 1*mm;
  EMModuleSizeZ    = 218*mm;
  HADModuleSizeZ   = 211*mm;
  SpacerThickness  = 3.5*mm;
  SpacerSizeY      = 21*mm;


  AbsorberThickness = 10.*mm;
  AbMedShieldThickness = 8.*mm;
  AbExShieldThickness = 1.*mm;
  GapThickness      = 3.5*mm;
  NbOfEMLayers      = 16;
  NbOfHADLayers     = 65;
  EMModuleSizeZ     = 218*mm;
  HADModuleSizeZ     = 211*mm;

  ModuleSizeY       = 700*mm;
  NbOfModules       = 3;
}


MyDetectorConstruction::~MyDetectorConstruction()
{}


G4VPhysicalVolume* MyDetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructCalorimeter();
}


void MyDetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String name, symbol;             //a=mass of a mole;
G4double a, z, density;            //z=mean number of protons;  
G4int iz, n;                       //iz=number of protons  in an isotope; 
                                   // n=number of nucleons in an isotope;

G4int ncomponents, natoms;
G4double abundance, fractionmass;
G4double temperature, pressure;

//
// define Elements
//

a = 1.0079*g/mole;
G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

a = 12.01*g/mole;
G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

a = 14.007*g/mole;
G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

a = 15.999*g/mole;
G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

a = 39.948*g/mole;
G4Element* Ar  = new G4Element(name="Argon"  ,symbol="Ar" , z= 18., a);

a = 207.19*g/mole;
G4Element* Pb  = new G4Element(name="Lead"  ,symbol="Pb" , z= 82., a);

a = 121.75*g/mole;
G4Element* Sb  = new G4Element(name="Antimony"  ,symbol="Sb" , z= 51., a);


//
// define simple materials
//

density     = universe_mean_density;    //from PhysicalConstants.h
pressure    = 3.e-18*pascal;
temperature = 2.73*kelvin;
G4Material* vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,
                   kStateGas,temperature,pressure);

density = 2.700*g/cm3;
a = 26.982*g/mole;
G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

// il numero di atomi e` stato normalizzato al numero di atomi della molecola!
// posso usare la formula chimica minimale, quindi moltiplicare per un 
// fattore comune fino ad ottenere interi.

density = 1.170*g/cm3;
G4Material* WLSMat = new G4Material(name="PMMAWLS", density, ncomponents=3);
WLSMat->AddElement(H, natoms=100);
WLSMat->AddElement(C, natoms=57);
WLSMat->AddElement(O, natoms=57);

//N.B. le frazioni di massa sono state normalizzate ad 1 rispetto all`input di fluka!
density = 11.3*g/cm3;
G4Material* ExLead = new G4Material(name="LeadSB", density, ncomponents=2);
ExLead->AddElement(Pb, fractionmass=0.96);
ExLead->AddElement(Sb, fractionmass=0.04);

density = 11.3*g/cm3;
G4Material* MedLead = new G4Material(name="LeadSB0", density, ncomponents=2);  
MedLead->AddElement(Pb, fractionmass=0.96);
MedLead->AddElement(Sb, fractionmass=0.04);

density = 0.001225*g/cm3;
G4Material* Air = new G4Material(name="Air", density, ncomponents=3);
Air->AddElement(N, fractionmass=0.7555795);
Air->AddElement(O, fractionmass=0.23158806);
Air->AddElement(Ar, fractionmass=0.012832444);

density = 1.044*g/cm3;
G4Material* SciMat = new G4Material(name="SCSN38", density, ncomponents=2);
SciMat->AddElement(H, natoms=1);
SciMat->AddElement(C, natoms=1);


//G4cout << *(G4Material::GetMaterialTable()) << endl;

  //default materials of the calorimeter
  GapMaterial      = Air;
  defaultMaterial  = Air;
  SpacerMaterial   = Al;
  RodMaterial      = SciMat;
  ScMaterial       = SciMat;

  WLSMaterial          = WLSMat;
  frontPlateMaterial   = Al;
  MedAbsMaterial      = MedLead;
  ExAbsMaterial       = ExLead;
}



  
G4VPhysicalVolume* MyDetectorConstruction::ConstructCalorimeter()
{
  // complete the Calor parameters definition 
  ComputeCalorParameters();
   
  //     
  // World
  //
  solidWorld = new G4Box("World",				//its name
                   WorldSizeX/2,WorldSizeY/2,WorldSizeZ/2);	//its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   defaultMaterial,	//vacuum
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
  
  //                               
  // Calorimeter
  //
      solidCalor = new G4Box("Calorimeter",		//its name
    		       ModuleSizeX/2,ModuleSizeY/2,CalorSizeZ/2);//size
    			     
      logicCalor = new G4LogicalVolume(solidCalor,	//its solid
      				       defaultMaterial,	//its material
      				       "Calorimeter");	//its name
    				       
      physiCalor = new G4PVPlacement(0,			//no rotation
      G4ThreeVector(ModuleSizeX/2-AlSizeX-SpaSizeX,0.,0.),
                                     "Calorimeter",	//its name
                                     logicCalor,	//its logical volume
                                     physiWorld,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number
  //                               
  // Module
  //
      solidModule = new G4Box("Module",		//its name
    		       ModuleSizeX/2,ModuleSizeY/2,ModuleSizeZ/2);//size
    			     
      logicModule = new G4LogicalVolume(solidModule,	        //its solid
				       defaultMaterial,      	//its material
      				       "Module");	        //its name


  for(int iMod=0;iMod<NbOfModules;iMod++)
    {
     new G4PVPlacement(0,
		       G4ThreeVector(0.,0.,
		       -ModuleSizeZ/2*(NbOfModules-1)+iMod*ModuleSizeZ),
		       "Module",
                       logicModule,
		       physiCalor,
		       false,
		       iMod);
    }


  //
  // EMWLS
  //
      solidEMWLS = new G4Box("EMWLS",		//its name
    		       EMWLSSizeX/2,ModuleSizeY/2,EMWLSSizeZ/2);//size
    			     
      logicEMWLS = new G4LogicalVolume(solidEMWLS,	        //its solid
				       WLSMaterial,           	//its material
      				       "EMWLS");	        //its name
      //right and left side of module
      for(int t=0; t<2; t++)
	{
	  int b=2*t-1;
	  new G4PVPlacement(0,			//no rotation
          G4ThreeVector(-ModuleSizeX/2+EMWLSSizeX/2+AlSizeX+SpaSizeX,
			0.,b*(-ModuleSizeZ/2+EMWLSSizeZ/2)),	//
                                     logicEMWLS,	//its logical volume 
                                     "EMWLS",       	//its name
                                     logicModule,	//its mother  volume
                                     false,		//no boolean operation
                                     t);	        //copy number
        }
                    						    
  //
  // HADWLS
  //
      solidHADWLS = new G4Box("HADWLS",		//its name
    		       HADWLSSizeX/2,ModuleSizeY/2,HADWLSSizeZ/2);//size
    			     
      logicHADWLS = new G4LogicalVolume(solidHADWLS,	        //its solid
				       WLSMaterial,          	//its material
      				       "HADWLS");	        //its name
      //left and right side of module
      for(int u=0; u<2; u++)
	{
	  int c=2*u-1;
	  new G4PVPlacement(0,			//no rotation
          G4ThreeVector(ModuleSizeX/2-HADWLSSizeX/2,
			0.,c*(ModuleSizeZ/2-EMWLSSizeZ-HADWLSSizeZ/2)),
                                     logicHADWLS,	//its logical volume
                                     "HADWLS",       	//its name
                                     logicModule,	//its mother  volume
                                     false,		//no boolean operation
                                     u);		//copy number
	}		
   							   
  
  //
  // Al front plate
  //
      solidAl = new G4Box("AlFrontPlate",		//its name
    		       AlSizeX/2,ModuleSizeY/2,ModuleSizeZ/2);//size
    			     
      logicAl = new G4LogicalVolume(solidAl,	        //its solid
				    frontPlateMaterial,             	//its material
      				    "AlFrontPlate");	        //its name
    				       
      physiAl = new G4PVPlacement(0,			//no rotation
                G4ThreeVector(-ModuleSizeX/2+AlSizeX/2,0.,0.),	
                                     logicAl,	        //its logical volume
                                     "AlFrontPlate",    	//its name
                                     logicModule,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number
					
  //
  // front space with vacuum 
  //
      solidSpace = new G4Box("VacFrontSpace",		//its name
    		       SpaSizeX/2,ModuleSizeY/2,ModuleSizeZ/2);//size
    			     
      logicSpace = new G4LogicalVolume(solidSpace,	        //its solid
				       defaultMaterial,            //its material
      				    "VacFrontSpace");	        //its name
    				       
      physiSpace = new G4PVPlacement(0,			//no rotation
                G4ThreeVector(-ModuleSizeX/2+AlSizeX+SpaSizeX/2,0.,0.),
                                     logicSpace,	        //its logical vo
                                     "VacFrontSpace",    	//its name
                                     logicModule,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number
  




  // *******************************                              
  // ********    EMModule    *******
  // *******************************

      EMsolidModule = new G4Box("EMModule",		//its name
    		       EMModuleThickness/2,ModuleSizeY/2,EMModuleSizeZ/2);//size
    			     
      EMlogicModule = new G4LogicalVolume(EMsolidModule,	//its solid
				       defaultMaterial,      	//its material
      				       "EMModule");	        //its name
    				       
      EMphysiModule = new G4PVPlacement(0,     		//no rotation
      G4ThreeVector(-ModuleSizeX/2+AlSizeX+SpaSizeX+EMModuleThickness/2,
		    0.,0.), 
                                     EMlogicModule,	//its logical volume
                                     "EMModule",	//its name
                                     logicModule,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number
  //                                 
  // EMLayer
  //
      EMsolidLayer = new G4Box("EMLayer",			//its name
                       LayerThickness/2,ModuleSizeY/2,EMModuleSizeZ/2); //size
                       
      EMlogicLayer = new G4LogicalVolume(EMsolidLayer,	        //its solid
                                       defaultMaterial,	        //its material
                                       "EMLayer");	        //its name
      
      G4VPVParameterisation * EMLayersParam = 
		     new T36LayerParam(LayerThickness,NbOfEMLayers);

      EMphysiLayer = new G4PVParameterised("EMLayer",     //its name
      				   EMlogicLayer,	  //its logical volume
      				   EMphysiModule,	  //its mother
                                   kXAxis,		  //axis of param
                                   NbOfEMLayers,	  //number of param
				   EMLayersParam);        //param

  //                               
  // Pb EMAbsorber
  //
      EMsolidAbsorber = new G4Box("EMAbsorber",		//its name
                          AbsorberThickness/2,ModuleSizeY/2,EMModuleSizeZ/2); 
                          
      EMlogicAbsorber = new G4LogicalVolume(EMsolidAbsorber,    //its solid
      			                  defaultMaterial,     //its material
      			                  "EMAbsorber");        //its name
      			                  
      EMphysiAbsorber = new G4PVPlacement(0,		     //no rotation
      		    G4ThreeVector(-GapThickness/2,0.,0.),    //its position
                                        "EMAbsorber",        //its name
                                        EMlogicAbsorber,     //its logical volume
                                        EMphysiLayer,        //its mother
                                        false,               //no boulean operat
                                        0);                  //copy number

  //                               
  // Pb EMAbsorber Medium Shield
  //
      EMsolidAbMedShield = new G4Box("EMAbMedShield",		//its name
                          AbMedShieldThickness/2,ModuleSizeY/2,EMModuleSizeZ/2); 
                          
      EMlogicAbMedShield = new G4LogicalVolume(EMsolidAbMedShield,    //its solid
      			                  MedAbsMaterial,                    //its material
      			                  "EMAbMedShield");           //its name
      			                  
      EMphysiAbMedShield = new G4PVPlacement(0,		     //no rotation
					G4ThreeVector(0.,0.,0.),    //its position
                                        "EMAbMedShield",        //its name
                                        EMlogicAbMedShield,     //its logical volume
                                        EMphysiAbsorber,        //its mother
                                        false,               //no boulean operat
                                        0);                  //copy number

  //                               
  // 2 Pb EMAbsorber External Shields
  //
      EMsolidAbExShield = new G4Box("EMAbExShield",		//its name
                          AbExShieldThickness/2,ModuleSizeY/2,EMModuleSizeZ/2); 
                          
      EMlogicAbExShield = new G4LogicalVolume(EMsolidAbExShield,    //its solid
      			                  ExAbsMaterial,                   //its material
      			                  "EMAbExShield");          //its name
								      
      			                 
      for(int s=0; s<2; s++)
	{
         EMphysiAbExShield = new G4PVPlacement(0,		     //no rotation
	 G4ThreeVector((2*s-1)*(AbsorberThickness/2-AbExShieldThickness/2),0.,0.),    //its position
                                        "EMAbExShield",        //its name
                                        EMlogicAbExShield,     //its logical volume
                                        EMphysiAbsorber,        //its mother
                                        false,               //no boulean operat
                                        s);                  //copy number
								 
        }

  
  //                                 
  // EMGap
  //
      EMsolidGap = new G4Box("EMGap",
    			   GapThickness/2,ModuleSizeY/2,EMModuleSizeZ/2);
    			   
      EMlogicGap = new G4LogicalVolume(EMsolidGap,
      				     GapMaterial,
      				     "EMGap");
      				     
      EMphysiGap = new G4PVPlacement(0,                      //no rotation
               G4ThreeVector(AbsorberThickness/2,0.,0.),     //its position
                                   "EMGap",                  //its name
                                   EMlogicGap,               //its logical volume
                                   EMphysiLayer,             //its mother
                                   false,                    //no boulean operat
                                   0);                       //copy number
   

  //
  // 3 EMScPlates
  //
      EMsolidSc = new G4Box("EMSc",
    			   ScThickness/2,ScSizeY/2,EMModuleSizeZ/2);
    			   
      EMlogicSc = new G4LogicalVolume(EMsolidSc,
      				     ScMaterial,
      				     "EMSc");
      				     

      for(int i=0; i<+3; i++)
	{
	  EMphysiSc = new G4PVPlacement(0,                   //no rotation
	    G4ThreeVector(0.,(i-1)*(ScSizeY+2*RodRad),0.),   //its position
                                   "EMSc",                   //its name
                                   EMlogicSc,                //its logical volume
                                   EMphysiGap,               //its mother
                                   false,                    //no boulean operat
                                   i);
	}  


  //
  // EMRod						
  //
      EMsolidRod = new G4Tubs("EMRod",
    			   0.,RodRad,EMModuleSizeZ/2,0*deg,360*deg);
    			   
      EMlogicRod = new G4LogicalVolume(EMsolidRod,
      				     RodMaterial,
      				     "EMRod");
      				     
      for(int h=0; h<2; h++)
	{
	  new G4PVPlacement(0,                                 //no rotation
          G4ThreeVector(0.,(2*h-1)*(ScSizeY/2+RodRad),0.),     //its position
                                   "EMRod",                    //its name
                                   EMlogicRod,                 //its logical volume
                                   EMphysiGap,                 //its mother
                                   false,                      //no boulean operat
                                   h);   	
	}

  //
  // EMSpacer
  //
      EMsolidSpacer = new G4Box("EMSpacer",
    			   SpacerThickness/2,SpacerSizeY/2,EMModuleSizeZ/2);
    			   
      EMlogicSpacer = new G4LogicalVolume(EMsolidSpacer,
      				     SpacerMaterial,
      				     "EMSpacer");
      		
      for(int k=0; k<2; k++)
	{
          int v=2*k-1;
	  new G4PVPlacement(0,                       //no rotation
	  G4ThreeVector(0.,v*(3*ScSizeY/2+2*RodRad+SpacerSizeY/2),0.),   //its position
                                   "EMSpacer",                   //its name
                                   EMlogicSpacer,                //its logical volume
                                   EMphysiGap,                   //its mother
                                   false,                        //no boulean operat
                                   k);
	}



  // *******************************                              
  // ********    HADModule    *******
  // *******************************

      HADsolidModule = new G4Box("HADModule",		//its name
    		       HADModuleThickness/2,ModuleSizeY/2,HADModuleSizeZ/2);//size
    			     
      HADlogicModule = new G4LogicalVolume(HADsolidModule,	//its solid
				       defaultMaterial,      	//its material
      				       "HADModule");	        //its name
    				       
      HADphysiModule = new G4PVPlacement(0,			//no rotation
      G4ThreeVector(ModuleSizeX/2-HADModuleThickness/2-
		    AbsorberThickness,0.,0.),	//at (0,0,0)
                                     HADlogicModule,	//its logical volume
                                     "HADModule",	//its name
                                     logicModule,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number
  
  //                                 
  // HADLayer
  //
      HADsolidLayer = new G4Box("HADLayer",			//its name
                       LayerThickness/2,ModuleSizeY/2,HADModuleSizeZ/2); //size
                       
      HADlogicLayer = new G4LogicalVolume(HADsolidLayer,	        //its solid
                                       defaultMaterial,	        //its material
                                       "HADLayer");	        //its name
                                      
      G4VPVParameterisation * HADLayersParam = 
		     new T36LayerParam(LayerThickness,NbOfHADLayers);

      HADphysiLayer = new G4PVParameterised("HADLayer",     //its name
      				   HADlogicLayer,	  //its logical volume
      				   HADphysiModule,	  //its mother
                                   kXAxis,		  //axis of param
                                   NbOfHADLayers,	  //number of param
				   HADLayersParam);        //param
  //                               
  // Pb HADAbsorber
  //
      HADsolidAbsorber = new G4Box("HADAbsorber",		//its name
                          AbsorberThickness/2,ModuleSizeY/2,HADModuleSizeZ/2); 
                          
      HADlogicAbsorber = new G4LogicalVolume(HADsolidAbsorber,    //its solid
      			                  defaultMaterial,     //its material
      			                  "HADAbsorber");        //its name
      			                  
      HADphysiAbsorber = new G4PVPlacement(0,		     //no rotation
      		    G4ThreeVector(-GapThickness/2,0.,0.),    //its position
                                        "HADAbsorber",        //its name
                                        HADlogicAbsorber,     //its logical volume
                                        HADphysiLayer,        //its mother
                                        false,               //no boulean operat
                                        0);                  //copy number


  //                               
  // Pb HADAbsorber Medium Shield
  //
      HADsolidAbMedShield = new G4Box("HADAbMedShield",		//its name
                          AbMedShieldThickness/2,ModuleSizeY/2,HADModuleSizeZ/2); 
                          
      HADlogicAbMedShield = new G4LogicalVolume(HADsolidAbMedShield,    //its solid
      			                  MedAbsMaterial,                      //its material
      			                  "HADAbMedShield");            //its name
      			                  
      HADphysiAbMedShield = new G4PVPlacement(0,		     //no rotation
					G4ThreeVector(0.,0.,0.),    //its position
                                        "HADAbMedShield",        //its name
                                        HADlogicAbMedShield,     //its logical volume
                                        HADphysiAbsorber,        //its mother
                                        false,               //no boulean operat
                                        0);                  //copy number

  //                               
  // 2 Pb HADAbsorber External Shields
  //
      HADsolidAbExShield = new G4Box("HADAbExShield",		//its name
                          AbExShieldThickness/2,ModuleSizeY/2,HADModuleSizeZ/2); 
                          
      HADlogicAbExShield = new G4LogicalVolume(HADsolidAbExShield,    //its solid
      			                  ExAbsMaterial,             //its material
      			                  "HADAbExShield");           //its name
      			                 
      for(int hs=0; hs<2; hs++)
	{
         HADphysiAbExShield = new G4PVPlacement(0,		     //no rotation
	 G4ThreeVector((2*hs-1)*(AbsorberThickness/2-AbExShieldThickness/2),0.,0.),    //its position
                                        "HADAbExShield",        //its name
                                        HADlogicAbExShield,     //its logical volume
                                        HADphysiAbsorber,        //its mother
                                        false,               //no boulean operat
                                        hs);                  //copy number
								 
        }

  //                               
  // last Pb HADAbsorber
  //


       lastHADsolidAbsorber = new G4Box("lastHADAbsorber",             //its name
                           AbsorberThickness/2,ModuleSizeY/2,HADModuleSizeZ/2);

       lastHADlogicAbsorber = new G4LogicalVolume(lastHADsolidAbsorber,    //its solid
                           defaultMaterial,                                //its material
                           "lastHADAbsorber");                             //its name

       lastHADphysiAbsorber = new G4PVPlacement(0,		     //no rotation
      		    G4ThreeVector(ModuleSizeX/2-AbsorberThickness/2,0.,0.),    //its position
                                        lastHADlogicAbsorber,     //its logical volume
                                        "lastHADAbsorber",        //its name
                                        logicModule,        //its mother
                                        false,               //no boulean operat
                                        0);                  //copy number
              
  //                               
  // last Pb HADAbsorber Medium Shield
  //
      lastHADphysiAbMedShield = new G4PVPlacement(0,		     //no rotation
					G4ThreeVector(0.,0.,0.),    //its position
                                        "lastHADAbMedShield",        //its name
                                        HADlogicAbMedShield,     //its logical volume
                                        lastHADphysiAbsorber,        //its mother
                                        false,               //no boulean operat
                                        0);                  //copy number

  //                               
  // 2 last Pb HADAbsorber External Shields
  //
      for(int ls=0; ls<2; ls++)
	{
         lastHADphysiAbMedShield = new G4PVPlacement(0,		     //no rotation
	 G4ThreeVector((2*ls-1)*(AbsorberThickness/2-AbExShieldThickness/2),0.,0.),    //its position
                                        "lastHADAbExShield",        //its name
                                        HADlogicAbExShield,     //its logical volume
                                        lastHADphysiAbsorber,        //its mother
                                        false,               //no boulean operat
                                        s);                  //copy number
								 
        }                                              
  

  //                                 
  // HADGap
  //
      HADsolidGap = new G4Box("HADGap",
    			   GapThickness/2,ModuleSizeY/2,HADModuleSizeZ/2);
    			   
      HADlogicGap = new G4LogicalVolume(HADsolidGap,
      				     GapMaterial,
      				     "HADGap");
      				     
      HADphysiGap = new G4PVPlacement(0,                      //no rotation
               G4ThreeVector(AbsorberThickness/2,0.,0.),     //its position
                                   "HADGap",                  //its name
                                   HADlogicGap,               //its logical volume
                                   HADphysiLayer,             //its mother
                                   false,                    //no boulean operat
                                   0);                       //copy number
   

  //
  // 3 HADScPlates
  //
      HADsolidSc = new G4Box("HADSc",
    			   ScThickness/2,ScSizeY/2,HADModuleSizeZ/2);
    			   
      HADlogicSc = new G4LogicalVolume(HADsolidSc,
      				     ScMaterial,
      				     "HADSc");
      				     

      for(int m=0; m<+3; m++)
	{
	  HADphysiSc = new G4PVPlacement(0,                  //no rotation
	    G4ThreeVector(0.,(m-1)*(ScSizeY+2*RodRad),0.),   //its position
                                   "HADSc",                   //its name
                                   HADlogicSc,                //its logical volume
                                   HADphysiGap,               //its mother
                                   false,                    //no boulean operat
                                   m);
	}  


  //
  // HADRod						
  //
      HADsolidRod = new G4Tubs("HADRod",
    			   0.,RodRad,HADModuleSizeZ/2,0*deg,360*deg);
    			   
      HADlogicRod = new G4LogicalVolume(HADsolidRod,
      				     RodMaterial,
      				     "HADRod");
      				     
      for(int q=0; q<2; q++)
	{
	  new G4PVPlacement(0,                                 //no rotation
          G4ThreeVector(0.,(2*q-1)*(ScSizeY/2+RodRad),0.),     //its position
                                   "HADRod",                    //its name
                                   HADlogicRod,                 //its logical volume
                                   HADphysiGap,                 //its mother
                                   false,                      //no boulean operat
                                   q);   	
	}

  //
  // HADSpacer
  //
      HADsolidSpacer = new G4Box("HADSpacer",
    			   SpacerThickness/2,SpacerSizeY/2,HADModuleSizeZ/2);
    			   
      HADlogicSpacer = new G4LogicalVolume(HADsolidSpacer,
      				     SpacerMaterial,
      				     "HADSpacer");
      		
      for(int z=0; z<2; z++)
	{
          int y=2*z-1;
	  new G4PVPlacement(0,                       //no rotation
	  G4ThreeVector(0.,y*(3*ScSizeY/2+2*RodRad+SpacerSizeY/2),0.),   //its position
                                   "HADSpacer",                   //its name
                                   HADlogicSpacer,                //its logical volume
                                   HADphysiGap,                   //its mother
                                   false,                        //no boulean operat
                                   z);
	}
      /*
  //
  // Sensitive EM and HAD Detectors: Absorber and Gap
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  SDman->SetVerboseLevel(1);

  //for N02: N.B. non basta dichiarare la "madre" volume sensibile, bisogna 
  // dichiarare anche ogni figlio!!!

  if(!EMcalorimeterSD)
  {
    G4String EMcaloSDname = "/mydet/EMcalorimeter";
    EMcalorimeterSD = new T36EMCalorimeterSD(EMcaloSDname,this);
    SDman->AddNewDetector( EMcalorimeterSD );
  }

  if (EMlogicAbsorber)
      EMlogicAbsorber->SetSensitiveDetector(EMcalorimeterSD);
  if (EMlogicAbExShield)
      EMlogicAbExShield->SetSensitiveDetector(EMcalorimeterSD);
  if (EMlogicAbMedShield)
      EMlogicAbMedShield->SetSensitiveDetector(EMcalorimeterSD);

  if (EMlogicGap)
      EMlogicGap->SetSensitiveDetector(EMcalorimeterSD);
  if (EMlogicSc)
      EMlogicSc->SetSensitiveDetector(EMcalorimeterSD);
  if (EMlogicRod)
       EMlogicRod->SetSensitiveDetector(EMcalorimeterSD);
  if (EMlogicSpacer)
       EMlogicSpacer->SetSensitiveDetector(EMcalorimeterSD);

  if(!HADcalorimeterSD)
  {
    G4String HADcaloSDname = "/mydet/HADcalorimeter";
    HADcalorimeterSD = new T36HADCalorimeterSD(HADcaloSDname,this);
    SDman->AddNewDetector( HADcalorimeterSD );
  }

  if (HADlogicAbsorber)
      HADlogicAbsorber->SetSensitiveDetector(HADcalorimeterSD);
  if (HADlogicAbExShield)
      HADlogicAbExShield->SetSensitiveDetector(HADcalorimeterSD);
  if (HADlogicAbMedShield)
      HADlogicAbMedShield->SetSensitiveDetector(HADcalorimeterSD);

  if (HADlogicGap)
      HADlogicGap->SetSensitiveDetector(HADcalorimeterSD);
  if (HADlogicSc)
      HADlogicSc->SetSensitiveDetector(HADcalorimeterSD);
  if (HADlogicRod)
       HADlogicRod->SetSensitiveDetector(HADcalorimeterSD);
  if (HADlogicSpacer)
       HADlogicSpacer->SetSensitiveDetector(HADcalorimeterSD);
      
  if (lastHADlogicAbsorber)
      lastHADlogicAbsorber->SetSensitiveDetector(HADcalorimeterSD);
      */

  /*

  //                                        
  // Visualization attributes
  //
								 
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  logicCalor->SetVisAttributes(simpleBoxVisAtt);
  
  */
  //
  //always return the physical World
  //
  return physiWorld;
}

void MyDetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n The calorimeter is made of " << NbOfEMLayers << "EM layers and "
       << NbOfHADLayers << "HAD layers of: [ "
       << AbsorberThickness/mm << "mm of " << MedAbsMaterial->GetName()
       << " + "
       << GapThickness/mm << "mm of " << GapMaterial->GetName() << " ] "
       << endl;
}

void MyDetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(magField) delete magField;		//delete the existing magn field
  
  if(fieldValue!=0.)			// create a new one if non nul
  { magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));        
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  } else {
    magField = NULL;
    fieldMgr->SetDetectorField(magField);
  }
}

  


