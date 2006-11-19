/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//////////////////////////////////////////////////////////////////////
//                                                                  //
//  (V-zero) detector  version 6  as designed by the Lyon group     //
//   All comments should be sent to Brigitte CHEYNIS :              //
//                                  b.cheynis@ipnl.in2p3.fr         // 
//   Geometry of september 2005 done with ROOT geometrical modeler  //                                  //
//   V0R (now V0C) sits between Z values  -89.5 and  -84.8 cm       //
//   V0L (now V0A) sits between Z values +339.0 and +341.0 cm       //
//   New coordinate system has been implemented in october 2003     //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "AliVZEROv6.h"
 
ClassImp(AliVZEROv6)

//_____________________________________________________________________________
AliVZEROv6:: AliVZEROv6():AliVZERO(),
   fCellId(0),
   fTrackPosition(),
   fTrackMomentum(), 
   fV0CHeight1(2.5), 
   fV0CHeight2(4.4), 
   fV0CHeight3(7.4), 
   fV0CHeight4(12.5),
   fV0CRMin(4.6), 
   fV0CRBox(38.0),
   fV0CLidThickness(0.30),
   fV0CCellThickness(2.00),
   fV0CBoxThickness(4.70),
   fV0COffsetFibers(1.0),
   fV0AHeight1(3.3), 
   fV0AHeight2(6.2), 
   fV0AHeight3(8.9), 
   fV0AHeight4(20.9),
   fV0ARMin(4.30),
   fV0ACellThickness(2.00),
   fLightYield(93.75),
   fLightAttenuation(0.05),
   fnMeters(15.0),
   fFibToPhot(0.3),
   fVersion(6)
{
// Standard default constructor 
}

//_____________________________________________________________________________
AliVZEROv6::AliVZEROv6(const char *name, const char *title):
   AliVZERO(name,title),
   fCellId(0),
   fTrackPosition(),
   fTrackMomentum(), 
   fV0CHeight1(2.5), 
   fV0CHeight2(4.4), 
   fV0CHeight3(7.4), 
   fV0CHeight4(12.5),
   fV0CRMin(4.6), 
   fV0CRBox(38.0),
   fV0CLidThickness(0.30),
   fV0CCellThickness(2.00),
   fV0CBoxThickness(4.70),
   fV0COffsetFibers(1.0),
   fV0AHeight1(3.3), 
   fV0AHeight2(6.2), 
   fV0AHeight3(8.9), 
   fV0AHeight4(20.9),
   fV0ARMin(4.30),
   fV0ACellThickness(2.00),
   fLightYield(93.75),
   fLightAttenuation(0.05),
   fnMeters(15.0),
   fFibToPhot(0.3),
   fVersion(6)
{

// Standard constructor for V-zero Detector  version 6

  AliDebug(2,"Create VZERO object ");
  
  fVersion            =     6;  // version number
  
// Parameters related to geometry :
// V0 part in front of muon arm absorber 

//   fV0CHeight1         =    2.5; // height of cell 1, in cm
//   fV0CHeight2         =    4.4; // height of cell 2, in cm
//   fV0CHeight3         =    7.4; // height of cell 3, in cm
//   fV0CHeight4         =   12.5; // height of cell 4, in cm
//   fV0CRMin            =    4.6; 
//   fV0CRBox            =   38.0; // outer radius of box, in cm
//   fV0CLidThickness    =   0.30; // thickness of Carbon lid
//   fV0CCellThickness   =   2.00; // thickness of elementary cell
//   fV0CBoxThickness    =   4.70; // thickness of V0C Box
//   fV0COffsetFibers    =    1.0; // offset to output fibers, in cm

// V0 part on the other side with respect to Interaction Point

//   fV0AHeight1         =    3.3; // height of cell 1, in cm
//   fV0AHeight2         =    6.2; // height of cell 2, in cm
//   fV0AHeight3         =    8.9; // height of cell 3, in cm
//   fV0AHeight4         =   20.9; // height of cell 4, in cm
//   fV0ARMin            =   4.30; 
//   fV0ACellThickness   =   2.00; // thickness of elementary cell  
//   
// Parameters related to light output :
         
//   fLightYield         =  93.75; // Light yield in BC408 (93.75 eV per photon)
//   fLightAttenuation   =   0.05; // Light attenuation in fiber (0.05 per meter)
//   fnMeters            =   15.0; // Number of meters of clear fibers to PM
//   fFibToPhot          =    0.3; // Attenuation at fiber-photocathode interface
}
     
//_____________________________________________________________________________

void AliVZEROv6::BuildGeometry()
{ 
          
}
            
//_____________________________________________________________________________
void AliVZEROv6::CreateGeometry()
{
  
// Constructs TGeo geometry 

  const int kColorVZERO  = kGreen;  
  
  AliDebug(2,"VZERO ConstructGeometry");
  
//  TGeoMedium  *medAir = gGeoManager->GetMedium("VZERO_Air"); 
  TGeoMedium  *medAlu = gGeoManager->GetMedium("VZERO_Aluminum");
  TGeoMedium  *medCar = gGeoManager->GetMedium("VZERO_Carbon");
  TGeoMedium  *medSci = gGeoManager->GetMedium("VZERO_Scintillator");
    
  TGeoVolume *top = gGeoManager->GetVolume("ALIC");
  
  Float_t  heightRight, r4Right;
  
  Float_t  zdet   =    90.0 - 0.5 - fV0CBoxThickness/2.0;
  heightRight     =    fV0CHeight1 + fV0CHeight2 + fV0CHeight3 + fV0CHeight4;
  r4Right         =    fV0CRMin + heightRight + 3.0*0.2;  // 3 spacings of 2mm between rings

// Creation of assembly V0RI - right part - :

  TGeoVolume *v0RI = new TGeoVolumeAssembly("V0RI");  
  TGeoTranslation *tr1 = new TGeoTranslation(0.,0.,-zdet);
  top->AddNode(v0RI,1,tr1);

// Creation of  carbon lids (3.0 mm thick) to keep V0C box shut :
    
  Float_t   partube[3];
  
  partube[0] =   fV0CRMin;
  partube[1] =   fV0CRBox;
  partube[2] =   fV0CLidThickness/2.0;
  
  TGeoTube   *sV0CA = new TGeoTube("V0CA", partube[0], partube[1], partube[2]);
  TGeoVolume *v0CA  = new TGeoVolume("V0CA",sV0CA,medCar);
  TGeoTranslation *tr2 = new TGeoTranslation(0.,0., fV0CBoxThickness/2.0-partube[2]);
  TGeoTranslation *tr3 = new TGeoTranslation(0.,0.,-fV0CBoxThickness/2.0+partube[2]);
  v0RI->AddNode(v0CA,1,tr2);
  v0RI->AddNode(v0CA,2,tr3);
  v0CA->SetLineColor(kYellow);
  
// Creation of aluminum rings 3.0 mm thick to maintain the v0RI pieces : 
 
  partube[0] =   fV0CRMin - 0.3;
  partube[1] =   fV0CRMin;
  partube[2] =   fV0CBoxThickness/2.0;
  
  TGeoTube   *sV0IR = new TGeoTube("V0IR", partube[0], partube[1], partube[2]);
  TGeoVolume *v0IR  = new TGeoVolume("V0IR",sV0IR,medAlu);
  v0RI->AddNode(v0IR,1,0);
  v0IR->SetLineColor(kYellow);
  
  partube[0] =   fV0CRBox;
  partube[1] =   fV0CRBox + 0.3; 
  partube[2] =   fV0CBoxThickness/2.0;

  TGeoTube   *sV0ER = new TGeoTube("V0ER", partube[0], partube[1], partube[2]);
  TGeoVolume *v0ER  = new TGeoVolume("V0ER",sV0ER,medAlu);
  v0RI->AddNode(v0ER,1,0);
  v0ER->SetLineColor(kYellow);
  
// Creation of assembly V0R0 of scintillator cells within one sector
 
  TGeoVolume *v0R0 = new TGeoVolumeAssembly("V0R0");  					  
  						 
// Elementary cell of ring 1  - right part - :
// (cells of ring 1 will be shifted by 2.0 cm backwards to output fibers)
  						  
  Float_t   r1Right =  fV0CRMin + fV0CHeight1;
  Float_t   offset  = fV0CBoxThickness/2.0 - fV0CLidThickness - fV0CCellThickness/2.0;   

  Float_t   partubs[5];   
     
  partubs[0]     =  fV0CRMin;
  partubs[1]     =  r1Right;
  partubs[2]     =  fV0CCellThickness/2.0;
  partubs[3]      =  90.0-22.5;
  partubs[4]      = 135.0-22.5;
  
  TGeoTubeSeg  *sV0R1 = new TGeoTubeSeg("V0R1", partubs[0], partubs[1], partubs[2], 
                                                partubs[3], partubs[4]);
  TGeoVolume   *v0R1  =  new TGeoVolume("V0R1",sV0R1,medSci);				       
  TGeoTranslation *tr4 = new TGeoTranslation(0.,0.,-offset);
  v0R0->AddNode(v0R1,1,tr4);
  v0R1->SetLineColor(kColorVZERO);

// Elementary cell of ring 2 - right part - :
// (cells of ring 2 will be shifted by 1.0 cm backwards to output fibers)

  Float_t   r2Right  =  r1Right + fV0CHeight2;  

  partubs[0]     =  r1Right;  //  must be equal to 7.1
  partubs[1]     =  r2Right;  //  must be equal to 11.5
  TGeoTubeSeg *sV0R2 = new TGeoTubeSeg("V0R2", partubs[0], partubs[1], partubs[2], 
                                               partubs[3], partubs[4]);
  TGeoVolume  *v0R2  = new TGeoVolume("V0R2",sV0R2,medSci);						      
  TGeoTranslation *tr5 = new TGeoTranslation(0.0,0.2,-offset + fV0COffsetFibers);					      
  v0R0->AddNode(v0R2,1,tr5);
  v0R2->SetLineColor(kColorVZERO);
   
// Ring 3 - right part -  :

//  Float_t   x = TMath::ATan(1.0/156.0) * ((180./TMath::Pi()));
  
  r2Right  =  r2Right + 0.2;
  Float_t   r3Right  =  r2Right + fV0CHeight3;     
//  printf(" r2 = %f, r3 = %f \n\n", r2Right,r3Right); 
  
  partubs[0]     =  r2Right;  //  must be equal to 11.7
  partubs[1]     =  r3Right;  //  must be equal to 19.1
  partubs[3]     =  90.0-22.5;
  partubs[4]     = 112.5-22.5;
  
  TGeoTubeSeg *sV0R3 = new TGeoTubeSeg("V0R3", partubs[0], partubs[1], partubs[2], 
                                               partubs[3], partubs[4]); 
  TGeoVolume  *v0R3  = new TGeoVolume("V0R3",sV0R3,medSci);					      
  TGeoTranslation *tr6 = new TGeoTranslation(0.,0.2,-offset + 2.0*fV0COffsetFibers);					       
  v0R0->AddNode(v0R3,1,tr6);
  v0R3->SetLineColor(kColorVZERO);
 
  partubs[3]     = 112.5-22.5;
  partubs[4]     = 135.0-22.5;
  
  TGeoTubeSeg *sV0R4 = new TGeoTubeSeg("V0R4", partubs[0], partubs[1], partubs[2], 
                                               partubs[3], partubs[4]);  
  TGeoVolume  *v0R4  = new TGeoVolume("V0R4",sV0R4,medSci);					       					      
  v0R0->AddNode(v0R4,1,tr6);
  v0R4->SetLineColor(kColorVZERO);
  
// Ring 4 - right part -  : 

  Float_t x = TMath::ATan(3.5/257.5) * ((180./TMath::Pi()));
  r3Right = r3Right + 0.2 + 0.2;   // + 0.2 because no shift in translation here !!
   
  partubs[0]     =  r3Right;  //  must be equal to 19.5
  partubs[1]     =  r4Right;  //  must be equal to 32.0
  partubs[3]     =  90.0-22.5+x;
  partubs[4]     = 112.5-22.5-x;
  
  TGeoTubeSeg *sV0R5 = new TGeoTubeSeg("V0R5", partubs[0], partubs[1], partubs[2], 
                                               partubs[3], partubs[4]);
  TGeoVolume  *v0R5  = new TGeoVolume("V0R5",sV0R5,medSci);
  TGeoTranslation *tr7 = new TGeoTranslation(0.,0.0,-offset + 2.0*fV0COffsetFibers);					      
  v0R0->AddNode(v0R5,1,tr7);
  v0R5->SetLineColor(kColorVZERO);
  
  partubs[3]     = 112.5-22.5+x;
  partubs[4]     = 135.0-22.5-x;
  
  TGeoTubeSeg *sV0R6 = new TGeoTubeSeg("V0R6", partubs[0], partubs[1], partubs[2], 
                                               partubs[3], partubs[4]);
  TGeoVolume  *v0R6  = new TGeoVolume("V0R6",sV0R6,medSci);						      
  v0R0->AddNode(v0R6,1,tr7);
  v0R6->SetLineColor(kColorVZERO);
  
  Float_t  phi;
  Float_t  phiDeg= 180./4.;
    
  Int_t    nsecR = 1;     // number of sectors in right part of V0
  Int_t    ncellsR;       // number of scintillating cells 
 
  for (phi = 22.5; phi < 360.0; phi = phi + phiDeg)
  
  {              
    TGeoRotation  *rot1 = new TGeoRotation("rot1", 90.0, +phi, 90., 90.+phi, 0.0, 0.0 ); 
    
    v0RI->AddNode(v0R0,nsecR,rot1);    
    nsecR++;        
  } 
     
  ncellsR = (nsecR - 1) * 6;    // 6 cells per sector (2 cells in  ring 3 and 4)  
  AliInfo(Form("Number of cells on Right side  - V0C =   %d",  ncellsR)); 
  
// Creation of assembly v0LE - left part - :
// Entrance face at  +339.0 cm  (new coordinate system) ...
          
  Float_t   heightLeft  = fV0AHeight1 + fV0AHeight2 + fV0AHeight3 + fV0AHeight4;   
  Float_t   r4Left      = fV0ARMin + heightLeft; 
 
  TGeoVolume *v0LE = new TGeoVolumeAssembly("V0LE"); 
   
  TGeoTranslation *tr8 = new TGeoTranslation(0.,0.,339.0 + fV0ACellThickness/2.0);
  top->AddNode(v0LE,1,tr8);
  
// Creation of assembly V0L0 of scintillator cells within one sector 
  
  TGeoVolume *v0L0 = new TGeoVolumeAssembly("V0L0");  					  
   
  Float_t   offsetLeft;
  offsetLeft    = - fV0ACellThickness/2.0; 

  Float_t   r1Left =  fV0ARMin + fV0AHeight1;        
      
  partubs[0]     =  fV0ARMin;
  partubs[1]     =  r1Left;
  partubs[2]      =  fV0ACellThickness/2.0;
  partubs[3]      =  90.0-22.5;
  partubs[4]      = 135.0-22.5;
  
  TGeoTubeSeg  *sV0L1 = new TGeoTubeSeg("V0L1", partubs[0], partubs[1], partubs[2], 
                                                partubs[3], partubs[4]);
  TGeoVolume   *v0L1  =  new TGeoVolume("V0L1",sV0L1,medSci);				       
  v0L0->AddNode(v0L1,1,gGeoIdentity);
  v0L1->SetLineColor(kColorVZERO);
  v0L1->SetVisibility(kTRUE);
  	 
  Float_t   r2Left =  r1Left + fV0AHeight2;       
  
  partubs[0]     =  r1Left;
  partubs[1]     =  r2Left;
  
  TGeoTubeSeg  *sV0L2 = new TGeoTubeSeg("V0L2", partubs[0], partubs[1], partubs[2], 
                                                partubs[3], partubs[4]);
  TGeoVolume   *v0L2  =  new TGeoVolume("V0L2",sV0L2,medSci);				       
  v0L0->AddNode(v0L2,1,gGeoIdentity);
  v0L2->SetLineColor(kColorVZERO);
  v0L2->SetVisibility(kTRUE);

  Float_t   r3Left =  r2Left + fV0AHeight3; 
   
  partubs[0]     =  r2Left;
  partubs[1]     =  r3Left;
  
  TGeoTubeSeg  *sV0L3 = new TGeoTubeSeg("V0L3", partubs[0], partubs[1], partubs[2], 
                                                partubs[3], partubs[4]);
  TGeoVolume   *v0L3  =  new TGeoVolume("V0L3",sV0L3,medSci);				       
  v0L0->AddNode(v0L3,1,gGeoIdentity);
  v0L3->SetLineColor(kColorVZERO);
  v0L3->SetVisibility(kTRUE);

  partubs[0]     =  r3Left;
  partubs[1]     =  r4Left;

  TGeoTubeSeg  *sV0L4 = new TGeoTubeSeg("V0L4", partubs[0], partubs[1], partubs[2], 
                                                partubs[3], partubs[4]);
  TGeoVolume   *v0L4  =  new TGeoVolume("V0L4",sV0L4,medSci);				       
  v0L0->AddNode(v0L4,1,gGeoIdentity);
  v0L4->SetLineColor(kColorVZERO);
  v0L4->SetVisibility(kTRUE);

  Int_t    nsecL = 1;     // number of sectors in left part of V0
  Int_t    ncellsL;       // number of scintillating cells 
  
  for (phi = 22.5; phi < 360.0; phi = phi + phiDeg)
  
  {                  
    TGeoRotation  *rot1 = new TGeoRotation("rot1", 90.0, +phi, 90., 90.+phi, 0.0, 0.0 ); 
    v0LE->AddNode(v0L0,nsecL,rot1);    
    nsecL++;        
  } 
     
  ncellsL = (nsecL - 1) * 4;    // 4 cells per sector
  AliInfo(Form("Number of cells on Left  side  - V0A =   %d\n",  ncellsL));

  gGeoManager->SetTopVolume(top); 
  gGeoManager->CloseGeometry();  
//  gGeoManager-> SetVisLevel(4);
}  
    
//_____________________________________________________________________________
void AliVZEROv6::CreateMaterials()
{

// Creates materials used for geometry 

   AliDebug(2,"Create materials");

//   Int_t  *idtmed = fIdtmed->GetArray()-2999;
      
   Int_t     fieldType       = gAlice->Field()->Integ();     // Field type 
   Double_t  maxField        = gAlice->Field()->Max();       // Field max.
   Double_t  maxBending      = 0;     // Max Angle
   Double_t  maxStepSize     = 0.001; // Max step size 
   Double_t  maxEnergyLoss   = 1;     // Max Delta E
   Double_t  precision       = 0.001; // Precision
   Double_t  minStepSize     = 0.001; // Minimum step size 
   Int_t     id;
   Double_t  a, z, density, radLength, absLength; 
   Float_t   tmaxfd, stemax, deemax, epsil, stmin;
   
   a = 0.0; z = 0.0; 
   density    = 0.0;
   radLength  = 0.0; 
   absLength  = 999.0;
   tmaxfd     = 10.;
   stemax     = 0.1;
   deemax     = 0.1;     
   epsil      = 0.001;
   stmin      = 0.001;
   
// Parameters  for Air (=  0.01% C + 75% N + 23% O + 1% Ar )

    Float_t aa[] = { 12.0107, 14.0067,   15.9994,  39.948 };
    Float_t za[] = {  6.,      7.,       8.,       18. };
    Float_t wa[] = { 0.000124, 0.755267, 0.231781, 0.012827 }; 
    density      = 0.00120479;
    maxBending   = 1;
    maxStepSize  = .001;
    precision    = .001;
    minStepSize  = .001;
    id           = 1;
    AliMixture(id, "Air", aa, za, density, 4, wa);
    AliMedium(id, "Air", id, 1, fieldType, maxField, maxBending,
		         maxStepSize, maxEnergyLoss, precision, minStepSize);
			
// Parameters  for Aluminum
 
    a = 26.98; 
    z = 13.00;
    density    = 2.7;
    radLength  = 8.9;
    maxBending  = 10;
    maxStepSize = .01;
    precision   = .003;
    minStepSize = .003;
    id = 2;
    AliMaterial( id, "Aluminum", a, z, density, radLength, 37.2, 0, 0);
    AliMedium(id, "Aluminum", id, 1, fieldType, maxField, maxBending,
		              maxStepSize, maxEnergyLoss, precision, minStepSize);
		    
// Parameters  for Carbon 

    a = 12.01; 
    z =  6.00;
    density     = 2.265;
    radLength   = 18.8;
    id = 3;
    AliMaterial(id, "Carbon",  a, z, density, radLength, 49.9, 0, 0);
    AliMedium(id,   "Carbon", id, 1, fieldType, maxField, maxBending,
		              maxStepSize, maxEnergyLoss, precision, minStepSize);
		    
// Parameters  for scintillator 

    Float_t as[] = { 1.00794, 12.011};
    Float_t zs[] = { 1.,  6.};
    Float_t ws[] = { 1.,  1.};
    density      = 1.032;
    maxBending   = 10;
    maxStepSize  = .01;
    precision    = .003;
    minStepSize  = .003;
    id           = 4;
    AliMixture(id, "Scintillator", as, zs, density, -2, ws);
    AliMedium(id,  "Scintillator", id, 1, fieldType, maxField, maxBending,
		                   maxStepSize,maxEnergyLoss,precision,minStepSize);

		                                             
}

//_____________________________________________________________________________
void AliVZEROv6::DrawModule() const
{

//  Drawing is done in DrawVZERO.C

   AliDebug(2,"DrawModule");
}


//_____________________________________________________________________________
void AliVZEROv6::DrawGeometry() 
{

//  Drawing of V0 geometry done in DrawV0.C

   AliDebug(2,"DrawGeometry");
 
//  Here is  DrawV0.C :

// void DrawV0()
// {
//    TGeoVolume *top = gGeoManager->GetMasterVolume();
//    gGeoManager->SetNsegments(80);
//    Int_t nd = top->GetNdaughters();
//    for (Int_t i=0; i<nd; i++) top->GetNode(i)->GetVolume()->InvisibleAll();
//    TGeoVolume *v0ri = gGeoManager->GetVolume("V0RI");  
//    TGeoVolume *v0le = gGeoManager->GetVolume("V0LE");
//    v0ri->SetVisibility(kTRUE);
//    v0ri->VisibleDaughters(kTRUE);
//    v0le->SetVisibility(kTRUE);
//    v0le->VisibleDaughters(kTRUE);
//    top->SetVisibility(kTRUE);
//    top->Draw();
// }
   
}

//_____________________________________________________________________________
void AliVZEROv6::Init()
{
// Initialises version of the VZERO Detector given in Config
// Just prints an information message
  
   AliInfo(Form("VZERO version %d initialized \n",IsVersion()));
   
   AliVZERO::Init();  
}

//_____________________________________________________________________________
void AliVZEROv6::StepManager()
{
 
// Step Manager, called at each step 
 
     Int_t     copy;
     static    Int_t   vol[4];
     static    Float_t hits[21];
     static    Float_t eloss, tlength;
     static    Int_t   nPhotonsInStep;
     static    Int_t   nPhotons; 
     static    Int_t   numStep;
     Float_t   ringNumber;
     Float_t   destep, step;
     
     numStep += 1; 
          
//   We keep only charged tracks :
     
     if ( !gMC->TrackCharge() || !gMC->IsTrackAlive() ) return; 

     vol[0]    = gMC->CurrentVolOffID(1, vol[1]);
     vol[2]    = gMC->CurrentVolID(copy);
     vol[3]    = copy;
     
     static Int_t idV0R1 = gMC->VolId("V0R1");
     static Int_t idV0L1 = gMC->VolId("V0L1");
     static Int_t idV0R2 = gMC->VolId("V0R2");
     static Int_t idV0L2 = gMC->VolId("V0L2");
     static Int_t idV0R3 = gMC->VolId("V0R3");
     static Int_t idV0L3 = gMC->VolId("V0L3");
     static Int_t idV0R4 = gMC->VolId("V0R4");
     static Int_t idV0L4 = gMC->VolId("V0L4");
     static Int_t idV0R5 = gMC->VolId("V0R5");
     static Int_t idV0R6 = gMC->VolId("V0R6");
   
     if      ( gMC->CurrentVolID(copy) == idV0R1 ||
               gMC->CurrentVolID(copy) == idV0L1 )
	       ringNumber = 1.0;
     else if ( gMC->CurrentVolID(copy) == idV0R2 ||
               gMC->CurrentVolID(copy) == idV0L2 ) 
	       ringNumber = 2.0;  
     else if ( gMC->CurrentVolID(copy) == idV0R3 ||
               gMC->CurrentVolID(copy) == idV0R4 ||
               gMC->CurrentVolID(copy) == idV0L3 )
	       ringNumber = 3.0;
     else if ( gMC->CurrentVolID(copy) == idV0R5 ||
               gMC->CurrentVolID(copy) == idV0R6 ||
               gMC->CurrentVolID(copy) == idV0L4 )
	       ringNumber = 4.0;	       
     else
     	       ringNumber = 0.0;
	       
 
     if  (  ringNumber > 0.5  ) { 
     
        destep    = gMC->Edep();
	step      = gMC->TrackStep();
        
	nPhotonsInStep  = Int_t(destep / (fLightYield *1e-9) );	
	nPhotonsInStep  = gRandom->Poisson(nPhotonsInStep);
	
	eloss    += destep;
	tlength  += step; 	 
	
        if  ( gMC->IsTrackEntering()  )  { 
	 
            nPhotons  =  nPhotonsInStep;       
	    gMC->TrackPosition(fTrackPosition);
	    gMC->TrackMomentum(fTrackMomentum);
	    
            Float_t pt  = TMath::Sqrt( fTrackMomentum.Px() * fTrackMomentum.Px() +
	                               fTrackMomentum.Py() * fTrackMomentum.Py() );
               
            hits[0]  = fTrackPosition.X();
            hits[1]  = fTrackPosition.Y();
            hits[2]  = fTrackPosition.Z();	 	 
	    hits[3]  = Float_t (gMC->TrackPid()); 

	    hits[4]  = gMC->TrackTime();
            hits[5]  = gMC->TrackCharge();
	    hits[6]  = fTrackMomentum.Theta()*TMath::RadToDeg();
	    hits[7]  = fTrackMomentum.Phi()*TMath::RadToDeg();
	    hits[8]  = ringNumber;
	 
	    hits[9]  = pt;
	    hits[10] = fTrackMomentum.P();
	    hits[11] = fTrackMomentum.Px();
	    hits[12] = fTrackMomentum.Py();
	    hits[13] = fTrackMomentum.Pz();
	    
	    TParticle *par = gAlice->GetMCApp()->Particle(gAlice->GetMCApp()->GetCurrentTrackNumber());
            hits[14] = par->Vx();
            hits[15] = par->Vy();
            hits[16] = par->Vz();
            
 	    tlength  = 0.0;
	    eloss    = 0.0;	    
         }
	 
	 nPhotons  = nPhotons + nPhotonsInStep;
	 
	 if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()){
	 
	 nPhotons  = nPhotons - Int_t((Float_t(nPhotons) * fLightAttenuation * fnMeters));	 
	 nPhotons  = nPhotons - Int_t( Float_t(nPhotons) * fFibToPhot);	 
	 
	 hits[17] =   eloss;
	 hits[18] = tlength;
	 hits[19] = nPhotons;
	 hits[20] = GetCellId (vol, hits); 
	  	 
         AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);
	 	 
	 tlength         = 0.0;
	 eloss           = 0.0; 
	 nPhotons        =   0;
	 nPhotonsInStep  =   0;
	 
	 numStep         =   0;  
	 } 
    }
      
}

//_____________________________________________________________________________
void AliVZEROv6::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  
//  Adds a VZERO hit
  
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliVZEROhit(fIshunt,track,vol,hits);
}

//_____________________________________________________________________________
void AliVZEROv6::AddDigits(Int_t *tracks, Int_t* digits) 
{

//  Adds a VZERO digit

   TClonesArray  &ldigits = *fDigits;
   new(ldigits[fNdigits++]) AliVZEROdigit(tracks, digits);
}

//_____________________________________________________________________________
void AliVZEROv6::MakeBranch(Option_t *option)
{
  
// Creates new branches in the current Root Tree
    
  char branchname[10];
  sprintf(branchname,"%s",GetName());
  AliDebug(2,Form("fBufferSize = %d",fBufferSize));
  
  const char *cH = strstr(option,"H");
  
  if (fHits   && TreeH() && cH) {
    TreeH()->Branch(branchname,&fHits, fBufferSize);
    AliDebug(2,Form("Making Branch %s for hits",branchname));
  }     

  const char *cD = strstr(option,"D");
  
  if (fDigits   && fLoader->TreeD() && cD) {
    fLoader->TreeD()->Branch(branchname,&fDigits, fBufferSize);
    AliDebug(2,Form("Making Branch %s for digits",branchname));
  }  
   
}

//_____________________________________________________________________________
Int_t AliVZEROv6::GetCellId(Int_t *vol, Float_t *hits) 
{

  //   Returns Id of scintillator cell
  //   Right side from  0 to 47
  //   Left  side from 48 to 95
  
  //   hits[8] = ring number (1 to 4)
  //   vol[1]  = copy number (1 to 8)

   Int_t index      = vol[1];
   Int_t ringNumber = Int_t(hits[8]);   
   fCellId          = 0;
   
//    cout << "volID = " << vol[0] << "  copy = " << vol[1] << endl;
//    cout << "X     = " << hits[0] << "    Y = " << hits[1] << endl;
   
   Float_t phi = Float_t(TMath::ATan2(Double_t(hits[1]),Double_t(hits[0])) ); 
   Float_t kRaddeg = 180.0/TMath::Pi();
   phi = kRaddeg * phi;
    
   if (index < 7) index = index + 8; 
   
   if (hits[2] < 0.0) { 
      if(ringNumber < 3) {
         index = (index - 7) + ( ( ringNumber - 1 ) * 8);}
      else if(ringNumber >= 3){ 
       if(gMC->CurrentVolID(vol[1]) == gMC->VolId("V0R3")|| 
          gMC->CurrentVolID(vol[1]) == gMC->VolId("V0R5") ) 
         {index = (index*2 - 14) + ( ( ringNumber - 2 ) * 16); }
       if(gMC->CurrentVolID(vol[1]) == gMC->VolId("V0R4")||
          gMC->CurrentVolID(vol[1]) == gMC->VolId("V0R6") ) 
         {index = (index*2 - 13) + ( ( ringNumber - 2 ) * 16); }
      }
      fCellId   = index;           
   }
           
   else if (hits[2] > 0.0){
      index = (index - 7 + 48) + ( ( ringNumber - 1 ) * 8);
      fCellId   = index;}
             
//    cout << " ring   = " << ringNumber << " phi = "<<  phi << endl; 
//    cout << " cellID = " << fCellId <<  endl;
//    cout <<  "**********" << endl;         
           
   return fCellId;
   
   
}
