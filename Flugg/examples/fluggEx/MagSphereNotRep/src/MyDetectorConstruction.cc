// $Id$
// Flugg tag $Name$

//
// Example with Sphere and Tub layers for testing Magnetic Field in FLUGG 
// Sara Vanini, 11/02/00. Not replicated volumes!
//


#include "MyDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4IntersectionSolid.hh"


MyDetectorConstruction::MyDetectorConstruction()
:solidWorld(NULL),logicWorld(NULL),physiWorld(NULL),
 solidDet(NULL),logicDet(NULL),physiDet(NULL),
 solidLayer(NULL),logicLayer(NULL),physiLayer1(NULL),
 physiLayer2(NULL),physiLayer3(NULL),physiLayer4(NULL),
 physiLayer5(NULL),physiLayer6(NULL),
 solidTub(NULL),logicTub(NULL),physiTub(NULL),
 solidSphere(NULL),logicSphere(NULL),physiSphere(NULL),
 TubMaterial(NULL),SphereMaterial(NULL),defaultMaterial(NULL),
 magField(NULL)
{
  // default parameter values of the calorimeter
  WorldSize         = 200.*cm;
  TubRad            = 3.*mm;
  SphereRad         = 5.*mm;
  NbOfLayers        = 6;
  DetSizeX          = 5.1*cm;
  DetSizeZ          = 8.*cm;
  NbOfTubs          = 13;
  NbOfSpheres       = 40;
}



MyDetectorConstruction::~MyDetectorConstruction()
{}



G4VPhysicalVolume* MyDetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructDetector();
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

a = 1.01*g/mole;
G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

a = 12.01*g/mole;
G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

a = 14.01*g/mole;
G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

a = 16.00*g/mole;
G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

a = 28.09*g/mole;
G4Element* Si = new G4Element(name="Silicon",symbol="Si" , z= 14., a);

a = 55.85*g/mole;
G4Element* Fe = new G4Element(name="Iron"    ,symbol="Fe", z=26., a);

//
// define an Element from isotopes, by relative abundance 
//

G4Isotope* U5 = new G4Isotope(name="U235", iz=92, n=235, a=235.01*g/mole);
G4Isotope* U8 = new G4Isotope(name="U238", iz=92, n=238, a=238.03*g/mole);

G4Element* U  = new G4Element(name="enriched Uranium", symbol="U", ncomponents=2);
U->AddIsotope(U5, abundance= 90.*perCent);
U->AddIsotope(U8, abundance= 10.*perCent);

//
// define simple materials
//

density = 2.700*g/cm3;
a = 26.98*g/mole;
G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

density = 1.390*g/cm3;
a = 39.95*g/mole;
G4Material* lAr = new G4Material(name="liquidArgon", z=18., a, density);

density = 8.960*g/cm3;
a = 63.55*g/mole;
G4Material* Cu = new G4Material(name="Copper"   , z=29., a, density);

density = 11.35*g/cm3;
a = 207.19*g/mole;
G4Material* Pb = new G4Material(name="Lead"     , z=82., a, density);

//
// define a material from elements.   case 1: chemical molecule
//
 
density = 1.000*g/cm3;
G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
H2O->AddElement(H, natoms=2);
H2O->AddElement(O, natoms=1);

density = 1.032*g/cm3;
G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
Sci->AddElement(C, natoms=9);
Sci->AddElement(H, natoms=10);

density = 2.200*g/cm3;
G4Material* SiO2 = new G4Material(name="quartz", density, ncomponents=2);
SiO2->AddElement(Si, natoms=1);
SiO2->AddElement(O , natoms=2);

//
// define a material from elements.   case 2: mixture by fractional mass
//

density = 1.290*mg/cm3;
G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
Air->AddElement(N, fractionmass=0.7);
Air->AddElement(O, fractionmass=0.3);

//
// define a material from elements and/or others materials (mixture of mixtures)
//

density = 0.200*g/cm3;
G4Material* Aerog = new G4Material(name="Aerogel", density, ncomponents=3);
Aerog->AddMaterial(SiO2, fractionmass=62.5*perCent);
Aerog->AddMaterial(H2O , fractionmass=37.4*perCent);
Aerog->AddElement (C   , fractionmass= 0.1*perCent);

//
// examples of gas in non STP conditions
//

density     = 27.*mg/cm3;
pressure    = 50.*atmosphere;
temperature = 325.*kelvin;
G4Material* CO2 = new G4Material(name="CarbonicGas", density, ncomponents=2,
                                     kStateGas,temperature,pressure);
CO2->AddElement(C, natoms=1);
CO2->AddElement(O, natoms=2);
 
density     = 0.3*mg/cm3;
pressure    = 2.*atmosphere;
temperature = 500.*kelvin;
G4Material* steam = new G4Material(name="WaterSteam", density, ncomponents=1,
                                      kStateGas,temperature,pressure);
steam->AddMaterial(H2O, fractionmass=1.);

//
// examples of vacuum
//

density     = universe_mean_density;    //from PhysicalConstants.h
pressure    = 3.e-18*pascal;
temperature = 2.73*kelvin;
new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,
                   kStateGas,temperature,pressure);

density     = 1.e-5*g/cm3;
pressure    = 2.e-2*bar;
temperature = STP_Temperature;         //from PhysicalConstants.h
G4Material* beam = new G4Material(name="Beam", density, ncomponents=1,
                                      kStateGas,temperature,pressure);
beam->AddMaterial(Air, fractionmass=1.);

//G4cout << *(G4Material::GetMaterialTable()) << endl;

  //default materials of the calorimeter
  TubMaterial = Al;
  SphereMaterial = Sci;
  LayerMaterial = Pb;
  defaultMaterial = Air;
}


  
G4VPhysicalVolume* MyDetectorConstruction::ConstructDetector()
{
  // complete the Detector parameters definition 
  ComputeDetectorParameters();
   
  //     
  // World
  //
  solidWorld = new G4Box("World",				//its name
                   WorldSize/2,WorldSize/2,WorldSize/2);	//its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                  defaultMaterial,            	//its material
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
  
  //                               
  // Detector
  //  
  solidDet=NULL; logicDet=NULL; physiDet=NULL;
  solidLayer=NULL; logicLayer=NULL; physiLayer1=NULL;
  physiLayer2=NULL; physiLayer3=NULL; physiLayer4=NULL;
  physiLayer5=NULL; physiLayer6=NULL;
  solidTubLayer=NULL; logicTubLayer=NULL; physiTubLayer=NULL;
  
      solidDet = new G4Box("Detector",		//its name
    		       DetSizeX/2,DetSizeY/2,DetSizeZ/2);//size
    			     
      logicDet = new G4LogicalVolume(solidDet,	//its solid
      				       defaultMaterial,	//its material
      				       "Detector");	//its name
    				       
      physiDet = new G4PVPlacement(0,			//no rotation
                                  G4ThreeVector(2.55*cm,-0.3*cm,4.0*cm),
				  "Detector",	//its name
                                  logicDet,	        //its logical volume
                                  physiWorld,	//its mother  volume
                                  false,		//no boolean operation
                                  0);		//copy number
  //                                 
  // Layer
  //
      solidLayer = new G4Box("Layer",			//its name
                       DetSizeX/2,LayerThickness/2,DetSizeZ/2); //size
                       
      logicLayer = new G4LogicalVolume(solidLayer,	//its solid
                                       LayerMaterial,	//its material
                                       "Layer");	//its name
      /*
      if (NbOfLayers > 1)                                      
        physiLayer = new G4PVReplica("Layer",		//its name
      				     logicLayer,	//its logical volume
      				     physiDet,	        //its mother
                                     kYAxis,		//axis of replication
                                     NbOfLayers,	//number of replica
                                     LayerThickness);	//witdth of replica
      else
        physiLayer = new G4PVPlacement(0,		//no rotation
                                     G4ThreeVector(),	//at (0,0,0)
                                     "Layer",		//its name
                                     logicLayer,	//its logical volume
                                     physiDet,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number     
      */

        physiLayer1 = new G4PVPlacement(0,		//no rotation
                                     G4ThreeVector(0,
				     -DetSizeY/2+LayerThickness/2,0),
                                     "Layer 1",		//its name
                                     logicLayer,	//its logical volume
                                     physiDet,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number     

        physiLayer2 = new G4PVPlacement(0,		//no rotation
                                     G4ThreeVector(0,
                                     -DetSizeY/2+3*LayerThickness/2,0),
                                     "Layer 2",		//its name
                                     logicLayer,	//its logical volume
                                     physiDet,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number     

        physiLayer3 = new G4PVPlacement(0,		//no rotation
                                     G4ThreeVector(0,
				     -LayerThickness/2,0),
                                     "Layer 3",		//its name
                                     logicLayer,	//its logical volume
                                     physiDet,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number     

        physiLayer4 = new G4PVPlacement(0,		//no rotation
                                     G4ThreeVector(0,
                                     +LayerThickness/2,0),
                                     "Layer 4",		//its name
                                     logicLayer,	//its logical volume
                                     physiDet,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number     

        physiLayer5 = new G4PVPlacement(0,		//no rotation
                                     G4ThreeVector(0,
				     +3*LayerThickness/2,0),
                                     "Layer 5",		//its name
                                     logicLayer,	//its logical volume
                                     physiDet,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number     

        physiLayer6 = new G4PVPlacement(0,		//no rotation
                                     G4ThreeVector(0,
				     +5*LayerThickness/2,0),
                                     "Layer 6",		//its name
                                     logicLayer,	//its logical volume
                                     physiDet,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number     






      /*
  //                                 
  // Tub-Layer
  //
      solidTubLayer = new G4Box("Tub-Layer",			//its name
                       TubRad*NbOfTubs,TubRad,DetSizeX/2); //size
                       
      logicTubLayer = new G4LogicalVolume(solidTubLayer,	//its solid
                                       LayerMaterial,	//its material
                                       "Tub Layer");	//its name

      G4RotationMatrix * rm = new G4RotationMatrix();
      G4double phi = 90*deg;
      rm->rotateY(phi);

      physiTubLayer = new G4PVPlacement(rm,		//rotation
                                     G4ThreeVector(0,SphereRad,-0.1*cm),
                                     logicTubLayer,	//its logical volume
				     "Tub Layer",		//its name
                                     logicLayer,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number     
      */
  //                               
  // Tubs
  //
  solidTub=NULL; logicTub=NULL; physiTub=NULL;solidBoxIntersTub=NULL;  
  solidTubSeg=NULL; logicTubSeg=NULL; physiTubSeg=NULL; 
  
  if (TubRad > 0.) 
    { solidTub = new G4Tubs("Tub",		//its name
			    0,
			    TubRad,
			    DetSizeX/2,
			    0.*deg,
			    360.*deg); 
                          
      logicTub = new G4LogicalVolume(solidTub,    //its solid
      			             TubMaterial, //its material
      			             "Tub");      //its name

      G4RotationMatrix * rm = new G4RotationMatrix();
      G4double phi = 90*deg;
      rm->rotateY(phi);

      for(int tubCopy=0; tubCopy<NbOfTubs; tubCopy++)
      {
        physiTub = new G4PVPlacement(rm,		//rotation
               G4ThreeVector(0,SphereRad,
		         -TubRad*NbOfTubs+TubRad*(1+2*tubCopy)-0.1*cm),
                         logicTub,	        //its logical volume
			 "Tub",	                //its name
                         logicLayer,	        //its mother  volume
                         false,		        //no boolean operation
                         tubCopy);		        //copy number

      }


      /*
      physiTub = new G4PVReplica("Tub",		//its name
      				 logicTub,	//its logical volume
      				 physiTubLayer,	        //its mother
                                 kXAxis,		//axis of replication
                                 NbOfTubs,	//number of replica
                                 2*TubRad);	//witdth of replica
      */
      
  //Tub segments
      solidBoxIntersTub = new G4Box("Tub segment",  //its name
		       TubRad,
		       TubRad,
		       DetSizeX/2);

      solidTubSeg = new G4IntersectionSolid("Tub segment",
					    solidTub,
					    solidBoxIntersTub,
					    0,
					    G4ThreeVector(-0.4*cm,0,0));

      logicTubSeg = new G4LogicalVolume(solidTubSeg,    //its solid
      			             TubMaterial, //its material
      			             "Tub segment");      //its name

      physiTubSeg = new G4PVPlacement(rm,		//rotation
                                   G4ThreeVector(0,SphereRad,4.1*cm),
				   logicTubSeg,	        //its logical volume   
                                   "Tub segment",	//its name
                                   logicLayer,	        //its mother  volume
                                   false,		//no boolean operation
                                   0);		        //copy number   
      
    }
  
  //                                 
  // Spheres
  //
  solidSphere=NULL; logicSphere=NULL; physiSphere=NULL; 
  
  if (SphereRad > 0.)
    { solidSphere = new G4Sphere("Sphere",
			      0*cm,SphereRad,
			      0,360*deg,
			      0,180*deg);
    		
      logicSphere = new G4LogicalVolume(solidSphere,
      				     SphereMaterial,
      				     "Sphere");
      
      G4double serie, element, Xposition, Yposition, Zposition;
      for(int copyNo=0;copyNo<40;copyNo++)
	{
	 serie = int(copyNo/9);
         element = copyNo - serie * 9 + 1;
         Yposition = -TubRad;

	 if (element<5)
	   {
	     Xposition = -2.55*cm + 2*SphereRad * element;
	     Zposition = -4.0*cm + SphereRad + sqrt(3)*SphereRad*2*serie;
	   }
	 else 
	   {
	     Xposition = -2.55*cm + SphereRad + (2*SphereRad)*(element-5);
	     Zposition = -4.0*cm + SphereRad + sqrt(3)*SphereRad*(2*serie+1);
	   }     

         physiSphere = new G4PVPlacement(0,		//rotation
                              G4ThreeVector(Xposition,Yposition,Zposition),
			      logicSphere,	        //its logical volume
			      "Sphere",	                //its name
                              logicLayer,	        //its mother  volume
                              false,     		//no boolean operation
                              copyNo);		        //copy number   
      
	}

  //Sphere segments
      solidBoxIntersSph = new G4Box("Sphere segment solid",  //its name
                       SphereRad,
		       SphereRad,
		       SphereRad);


      solidSphSeg1 = new G4IntersectionSolid("Sphere segment 1",
                                            solidSphere,
					    solidBoxIntersSph,
					    0,
					    G4ThreeVector(SphereRad,0,0));

      logicSphSeg1 = new G4LogicalVolume(solidSphSeg1,
      				     SphereMaterial,
      				     "Sphere segment 1");

      solidSphSeg2 = new G4IntersectionSolid("Sphere segment 2",
					     solidSphere,
  					     solidBoxIntersSph,
 					     0,
					     G4ThreeVector(-0.4*cm,0,0));

      logicSphSeg2 = new G4LogicalVolume(solidSphSeg2,
      				     SphereMaterial,
      				     "Sphere segment 2");

      solidSphSeg3 = new G4IntersectionSolid("Sphere segment 3",
					     solidSphere,
					     solidBoxIntersSph,
					     0,
					     G4ThreeVector(-0.9*cm,0,0));

      logicSphSeg3 = new G4LogicalVolume(solidSphSeg3,
      				     SphereMaterial,
      				     "Sphere segment 3");

      solidSphSeg4 = new G4IntersectionSolid("Sphere segment 4", 
					     solidSphere,
					     solidBoxIntersSph,
					     0,
					     G4ThreeVector(0,0,
                                             8.0*cm-(9*sqrt(3)+2)*SphereRad));

      logicSphSeg4 = new G4LogicalVolume(solidSphSeg4,
      				     SphereMaterial,
      				     "Sphere segment 4");


      for(int s1=0;s1<5;s1++)
      {
        physiSphSeg1 = new G4PVPlacement(0,		//no rotation
                         G4ThreeVector(-DetSizeX/2,-TubRad,
			 -4.0*cm+SphereRad+sqrt(3)*SphereRad*2*s1),
                         logicSphSeg1,	        //its logical volume
			 "Sphere segments 1",	        //its name
                         logicLayer,	        //its mother  volume
                         false,		        //no boolean operation
                         s1);		        //copy number

        physiSphSeg2 = new G4PVPlacement(0,		//no rotation
                         G4ThreeVector(2.45*cm,-TubRad,
			 -4.0*cm+SphereRad+sqrt(3)*SphereRad*2*s1),
                         logicSphSeg2,	        //its logical volume
			 "Sphere segments 2",	        //its name
                         logicLayer,	        //its mother  volume
                         false,		        //no boolean operation
                         s1);		        //copy number
       
        physiSphSeg4 = new G4PVPlacement(0,		//no rotation
                         G4ThreeVector(
			 -2.55*cm+SphereRad+2*SphereRad*s1,
			 -TubRad,
			 -4.0*cm+SphereRad+sqrt(3)*SphereRad*9),
                         logicSphSeg4,	        //its logical volume
                         "Sphere segments 4",	        //its name
                         logicLayer,	        //its mother  volume
                         false,		        //no boolean operation
                         s1);		        //copy number
      }
      
      for(int s2=0;s2<4;s2++)
      {
        physiSphSeg3 = new G4PVPlacement(0,		//no rotation
                         G4ThreeVector(2.95*cm,
			 -TubRad,
			 -4.0*cm+SphereRad+sqrt(3)*SphereRad*(2*s2+1)),
                         logicSphSeg3,	        //its logical volume
                         "Sphere segments 3",	        //its name
                         logicLayer,	        //its mother  volume
                         false,		        //no boolean operation
                         s2);		        //copy number
			 }
    }

  
  //                                        
  // Visualization attributes
  //
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  logicDet->SetVisAttributes(simpleBoxVisAtt);

  //
  //always return the physical World
  //

  return physiWorld;
}



/*
void MyDetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
  G4FieldManager* fieldMgr = ptrGeoInit->getFieldManager();
    
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

  */
