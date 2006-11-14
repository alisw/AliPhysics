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

///////////////////////////////////////////////////////////////////////
//                                                                   //
//  (V-zero) detector  version 7 as designed by the Lyon and         //
//   Mexico groups and Carlos Perez Lara from Pontificia Universidad //
//   Catolica del Peru    				             // 
//   All comments should be sent to Brigitte CHEYNIS:                //
//                     b.cheynis@ipnl.in2p3.fr                       // 
//   Geometry of April 2006 done with ROOT geometrical modeler       //
//   V0R (now V0C) sits between Z values  -89.5 and  -84.8 cm        //
//   V0L (now V0A) sits between Z values +338.5 and +342.5 cm        //
//   New coordinate system has been implemented in october 2003      //
//                                                                   //
///////////////////////////////////////////////////////////////////////

// --- Standard libraries ---
#include <Riostream.h>

// --- ROOT libraries ---
#include <TClonesArray.h>
#include <TMath.h>
#include <TVirtualMC.h>
#include <TParticle.h>

#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include "TGeoTube.h"
#include "TGeoArb8.h"
#include "TGeoCompositeShape.h"

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliVZEROLoader.h"
#include "AliVZEROdigit.h"
#include "AliVZEROhit.h"
#include "AliVZEROv7.h"
#include "AliLog.h"
 
ClassImp(AliVZEROv7)

//_____________________________________________________________________________
AliVZEROv7:: AliVZEROv7():AliVZERO()
{
// Standard default constructor 
}

//_____________________________________________________________________________
AliVZEROv7::AliVZEROv7(const char *name, const char *title):AliVZERO(name,title)
{
// Standard constructor for V-zero Detector  version 7

  AliDebug(2,"Create VZERO object ");
  fVersion            =     7;  // version number

  // V0C Parameters related to geometry: All in cm
  fV0CHeight1         =    2.5; // height of cell 1
  fV0CHeight2         =    4.4; // height of cell 2
  fV0CHeight3         =    7.4; // height of cell 3
  fV0CHeight4         =   12.5; // height of cell 4
  fV0CRMin            =    4.6; // inner radius of box
  fV0CRBox            =   38.0; // outer radius of box
  fV0CLidThickness    =   0.30; // thickness of Carbon lid
  fV0CCellThickness   =   2.00; // thickness of elementary cell
  fV0CBoxThickness    =   4.70; // thickness of V0C Box
  fV0COffsetFibers    =    1.0; // offset to output fibers
  // V0C Parameters related to light output
  fV0CLightYield         =  93.75; // Light yield in BC408 (93.75 eV per photon)
  fV0CLightAttenuation   =   0.05; // Light attenuation in fiber (0.05 per meter)
  fV0CnMeters            =   15.0; // Number of meters of clear fibers to PM
  fV0CFibToPhot          =    0.3; // Attenuation at fiber-photocathode interface

  // V0A Parameters related to geometry: All in cm
  fV0AR0     =  4.2;  // Radius of hole
  fV0AR1     =  7.6;  // Maximun radius of 1st cell
  fV0AR2     = 13.8; // Maximun radius of 2nd cell
  fV0AR3     = 22.7; // Maximun radius of 3rd cell
  fV0AR4     = 41.3; // Maximun radius of 4th cell
  fV0AR5     = 43.3; // Radius circunscrite to innermost octagon
  fV0AR6     = 68.0; // Radius circunscrite to outtermost octagon
  fV0ASciWd  =  2.5;  // Scintillator thickness 
  fV0APlaWd  =  0.5;  // Plates thinckness
  fV0APlaAl  = 0.06; // Plates AlMg3 thinckness
  fV0AOctWd  = 0.75; // Innermost octagon thickness
  fV0AOctH1  =  1.0;  // Height of innermost octagon
  fV0AOctH2  =  2.0;  // Height of outtermost octagon
  fV0AFibRd  =  0.1;  // Radius of Fiber
  fV0AFraWd  =  0.2;  // Support Frame thickness
  fV0APMBWd  = 24.6;  // Width of PM Box
  fV0APMBHt  = 22.0;  // Height of PM Box
  fV0APMBTh  =  7.1;  // Thickness of PM Box
  fV0APMBWdW =  0.3;  // Thickness of PM Box Side1 Wall
  fV0APMBHtW =  1.0;  // Thickness of PM Box Side2 Wall
  fV0APMBThW =  0.3;  // Thickness of PM Box Top Wall
  fV0APMBAng = 30.0;  // Angle between PM Box and Support
  fV0APMTR1  = 2.44;  // PMT Glass
  fV0APMTR2  = 2.54;  // PMT Glass
  fV0APMTR3  = 2.54;  // PMT Cover
  fV0APMTR4  = 2.70;  // PMT Cover
  fV0APMTH   = 10.0;  // PMT Height
  fV0APMTB   =  1.0;  // PMT Basis
  fV0APlaEx  =  4.4;  // Plates Extension height
  fV0ABasHt  =  2.0;  // Basis Height
  // V0A Parameters related to light output
  fV0ALightYield         =  93.75;      // Light yield in BC404
  fV0ALightAttenuation   =   0.05;      // Light attenuation in WLS fiber, per meter
  fV0AnMeters            = fV0AR6*0.01; // Tentative value, in meters
  fV0AFibToPhot          =    0.3;      // Attenuation at fiber-photocathode interface
}
//_____________________________________________________________________________

void AliVZEROv7::BuildGeometry()
{ 
}
            
//_____________________________________________________________________________
void AliVZEROv7::CreateGeometry()
{
// Constructs TGeo geometry 

  AliDebug(2,"VZERO ConstructGeometry");
  TGeoVolume *top = gGeoManager->GetVolume("ALIC");

  ///////////////////////////////////////////////////////////////////////////
  // Construct the geometry of V0C Detector. Brigitte CHEYNIS
  
    const int kColorVZERO  = kGreen;
    TGeoMedium *medV0CAlu = gGeoManager->GetMedium("VZERO_V0CAlu");
    TGeoMedium *medV0CCar = gGeoManager->GetMedium("VZERO_V0CCar");
    TGeoMedium *medV0CSci = gGeoManager->GetMedium("VZERO_V0CSci");
    TGeoVolume *v0RI = new TGeoVolumeAssembly("V0RI");
    Float_t heightRight, r4Right;
    Float_t zdet = 90.0 - 0.5 - fV0CBoxThickness/2.0;
    heightRight  = fV0CHeight1 + fV0CHeight2 + fV0CHeight3 + fV0CHeight4;
    r4Right      = fV0CRMin + heightRight + 3.0*0.2; // 3 spacings of 2mm between rings

    // Creation of  carbon lids (3.0 mm thick) to keep V0C box shut :
    Float_t   partube[3];
    partube[0] =   fV0CRMin;
    partube[1] =   fV0CRBox;
    partube[2] =   fV0CLidThickness/2.0;
    TGeoTube   *sV0CA = new TGeoTube("V0CA", partube[0], partube[1], partube[2]);
    TGeoVolume *v0CA  = new TGeoVolume("V0CA",sV0CA,medV0CCar);
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
    TGeoVolume *v0IR  = new TGeoVolume("V0IR",sV0IR,medV0CAlu);
    v0RI->AddNode(v0IR,1,0);
    v0IR->SetLineColor(kYellow);
    partube[0] =   fV0CRBox;
    partube[1] =   fV0CRBox + 0.3; 
    partube[2] =   fV0CBoxThickness/2.0;
    TGeoTube   *sV0ER = new TGeoTube("V0ER", partube[0], partube[1], partube[2]);
    TGeoVolume *v0ER  = new TGeoVolume("V0ER",sV0ER,medV0CAlu);
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
    partubs[3]     =  90.0-22.5;
    partubs[4]     = 135.0-22.5;
    TGeoTubeSeg *sV0R1 = new TGeoTubeSeg("V0R1", partubs[0], partubs[1], partubs[2],
					 partubs[3], partubs[4]);
    TGeoVolume  *v0R1  = new TGeoVolume("V0R1",sV0R1,medV0CSci);				       
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
    TGeoVolume  *v0R2  = new TGeoVolume("V0R2",sV0R2,medV0CSci);
    TGeoTranslation *tr5 = new TGeoTranslation(0.0,0.2,-offset + fV0COffsetFibers);
    v0R0->AddNode(v0R2,1,tr5);
    v0R2->SetLineColor(kColorVZERO);

    // Ring 3 - right part -  :
    r2Right  =  r2Right + 0.2;
    Float_t   r3Right  =  r2Right + fV0CHeight3;     
    partubs[0]     =  r2Right;  //  must be equal to 11.7
    partubs[1]     =  r3Right;  //  must be equal to 19.1
    partubs[3]     =  90.0-22.5;
    partubs[4]     = 112.5-22.5;
    TGeoTubeSeg *sV0R3 = new TGeoTubeSeg("V0R3", partubs[0], partubs[1], partubs[2],
					 partubs[3], partubs[4]);
    TGeoVolume  *v0R3  = new TGeoVolume("V0R3",sV0R3,medV0CSci);
    TGeoTranslation *tr6 = new TGeoTranslation(0.,0.2,-offset + 2.0*fV0COffsetFibers);
    v0R0->AddNode(v0R3,1,tr6);
    v0R3->SetLineColor(kColorVZERO);
    partubs[3]     = 112.5-22.5;
    partubs[4]     = 135.0-22.5;
    TGeoTubeSeg *sV0R4 = new TGeoTubeSeg("V0R4", partubs[0], partubs[1], partubs[2],
					 partubs[3], partubs[4]);
    TGeoVolume  *v0R4  = new TGeoVolume("V0R4",sV0R4,medV0CSci);
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
    TGeoVolume  *v0R5  = new TGeoVolume("V0R5",sV0R5,medV0CSci);
    TGeoTranslation *tr7 = new TGeoTranslation(0.,0.0,-offset + 2.0*fV0COffsetFibers);					      
    v0R0->AddNode(v0R5,1,tr7);
    v0R5->SetLineColor(kColorVZERO);
    partubs[3]     = 112.5-22.5+x;
    partubs[4]     = 135.0-22.5-x;
    TGeoTubeSeg *sV0R6 = new TGeoTubeSeg("V0R6", partubs[0], partubs[1], partubs[2],
					 partubs[3], partubs[4]);
    TGeoVolume  *v0R6  = new TGeoVolume("V0R6",sV0R6,medV0CSci);
    v0R0->AddNode(v0R6,1,tr7);
    v0R6->SetLineColor(kColorVZERO);
    Float_t  phi;
    Float_t  phiDeg= 180./4.;
    Int_t    nsecR = 1;     // number of sectors in right part of V0
    for (phi = 22.5; phi < 360.0; phi = phi + phiDeg) {
      TGeoRotation  *rot1 = new TGeoRotation("rot1", 90.0, +phi, 90., 90.+phi, 0.0, 0.0 ); 
      v0RI->AddNode(v0R0,nsecR,rot1);    
      nsecR++;        
    }

  ///////////////////////////////////////////////////////////////////////////
  // Construct the geometry of V0A Detector. Carlos PEREZ, PUCP

    const int kV0AColorSci   = 5;
    const int kV0AColorPlaIn = 3;
    const int kV0AColorPlaOu = 41;
    const int kV0AColorOct   = 7;
    const int kV0AColorFra   = 6;
    const int kV0AColorFib   = 11;
    const int kV0AColorPMG   = 1;
    const int kV0AColorPMA   = 2;
    const int kV0AColorBas   = 20;
    TGeoMedium *medV0ASci = gGeoManager->GetMedium("VZERO_V0ASci");
    TGeoMedium *medV0APlaIn = gGeoManager->GetMedium("VZERO_V0APlaIn");
    TGeoMedium *medV0APlaOu = gGeoManager->GetMedium("VZERO_V0APlaOu");
    TGeoMedium *medV0ASup = gGeoManager->GetMedium("VZERO_V0ALuc");
    TGeoMedium *medV0AFra = gGeoManager->GetMedium("VZERO_V0ALuc");
    TGeoMedium *medV0AFib = gGeoManager->GetMedium("VZERO_V0AFib");
    TGeoMedium *medV0APMGlass = gGeoManager->GetMedium("VZERO_V0APMG");
    TGeoMedium *medV0APMAlum = gGeoManager->GetMedium("VZERO_V0APMA");
    TGeoMedium *medV0ABas = gGeoManager->GetMedium("VZERO_V0ALuc");
    double pi = TMath::Pi();
    double sin225   = TMath::Sin(pi/8.);
    double cos225   = TMath::Cos(pi/8.);
    double ctg225   = cos225/sin225;
    double sin45    = TMath::Sin(pi/4.); // lucky: Sin45=Cos45
    double v0APts[16];

    ////////////////////////////
    /// Definition of one sector
    TGeoVolume *v0ASec = new TGeoVolumeAssembly("V0ASec");
    
    /// For boolean sustraction
    double preShape = 0.2;
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0-fV0AFraWd/2.-preShape;  v0APts[1+8*i] = -preShape;
      v0APts[2+8*i] = fV0AR0-fV0AFraWd/2.-preShape;  v0APts[3+8*i] = fV0AFraWd/2.;
      v0APts[4+8*i] = fV0AR4+fV0AFraWd/2.+preShape;  v0APts[5+8*i] = fV0AFraWd/2.;
      v0APts[6+8*i] = fV0AR4+fV0AFraWd/2.+preShape;  v0APts[7+8*i] = -preShape;
    }
    new TGeoArb8("sV0ACha1",fV0ASciWd/1.5,v0APts);
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0*sin45-preShape;
      v0APts[1+8*i] = (fV0AR0-fV0AFraWd)*sin45-preShape;
      v0APts[2+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45-preShape;
      v0APts[3+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45;
      v0APts[4+8*i] = (fV0AR4+fV0AFraWd/2.)*sin45+preShape;
      v0APts[5+8*i] = (fV0AR4+fV0AFraWd/2.)*sin45+2.*preShape;
      v0APts[6+8*i] = (fV0AR4+fV0AFraWd)*sin45+preShape;
      v0APts[7+8*i] = fV0AR4*sin45+preShape;
    }
    new TGeoArb8("sV0ACha2", fV0ASciWd/2.+2.*preShape, v0APts);
    new TGeoCompositeShape("sV0ACha","sV0ACha1+sV0ACha2");

    /// Frame
    TGeoVolume *v0AFra = new TGeoVolumeAssembly("V0AFra");
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[3+8*i] = fV0AFraWd/2.;
      v0APts[4+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[5+8*i] = fV0AFraWd/2.;
      v0APts[6+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[7+8*i] = 0.;
    }
    TGeoArb8 *sV0AFraB1 = new TGeoArb8("sV0AFraB1",fV0ASciWd/2.,v0APts);
    TGeoVolume *v0AFraB1 = new TGeoVolume("V0AFraB1",sV0AFraB1,medV0AFra);
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0*sin45;
      v0APts[1+8*i] = (fV0AR0-fV0AFraWd)*sin45;
      v0APts[2+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45;
      v0APts[3+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45;
      v0APts[4+8*i] = (fV0AR4+fV0AFraWd/2.)*sin45;
      v0APts[5+8*i] = (fV0AR4+fV0AFraWd/2.)*sin45;
      v0APts[6+8*i] = (fV0AR4+fV0AFraWd)*sin45;
      v0APts[7+8*i] = fV0AR4*sin45;
    }
    TGeoArb8 *sV0AFraB2 = new TGeoArb8("sV0AFraB2", fV0ASciWd/2., v0APts);
    TGeoVolume *v0AFraB2 = new TGeoVolume("V0AFraB2",sV0AFraB2,medV0AFra);
    v0AFraB1->SetLineColor(kV0AColorFra); v0AFraB2->SetLineColor(kV0AColorFra);
    v0AFra->AddNode(v0AFraB1,1);
    v0AFra->AddNode(v0AFraB2,1);  // Prefer 2 GeoObjects insted of 3 GeoMovements
    new TGeoTubeSeg( "sV0AFraR1b", fV0AR0-fV0AFraWd/2.,
		     fV0AR0+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR2b", fV0AR1-fV0AFraWd/2.,
		     fV0AR1+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR3b", fV0AR2-fV0AFraWd/2.,
		     fV0AR2+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR4b", fV0AR3-fV0AFraWd/2.,
		     fV0AR3+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR5b", fV0AR4-fV0AFraWd/2.,
		     fV0AR4+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    TGeoCompositeShape *sV0AFraR1 = new TGeoCompositeShape("sV0AFraR1","sV0AFraR1b-sV0ACha");
    TGeoCompositeShape *sV0AFraR2 = new TGeoCompositeShape("sV0AFraR2","sV0AFraR2b-sV0ACha");
    TGeoCompositeShape *sV0AFraR3 = new TGeoCompositeShape("sV0AFraR3","sV0AFraR3b-sV0ACha");
    TGeoCompositeShape *sV0AFraR4 = new TGeoCompositeShape("sV0AFraR4","sV0AFraR4b-sV0ACha");
    TGeoCompositeShape *sV0AFraR5 = new TGeoCompositeShape("sV0AFraR5","sV0AFraR5b-sV0ACha");
    TGeoVolume *v0AFraR1 = new TGeoVolume("V0AFraR1",sV0AFraR1,medV0AFra);
    TGeoVolume *v0AFraR2 = new TGeoVolume("V0AFraR2",sV0AFraR2,medV0AFra);
    TGeoVolume *v0AFraR3 = new TGeoVolume("V0AFraR3",sV0AFraR3,medV0AFra);
    TGeoVolume *v0AFraR4 = new TGeoVolume("V0AFraR4",sV0AFraR4,medV0AFra);
    TGeoVolume *v0AFraR5 = new TGeoVolume("V0AFraR5",sV0AFraR5,medV0AFra);
    v0AFraR1->SetLineColor(kV0AColorFra); v0AFraR2->SetLineColor(kV0AColorFra);
    v0AFraR3->SetLineColor(kV0AColorFra); v0AFraR4->SetLineColor(kV0AColorFra);
    v0AFraR5->SetLineColor(kV0AColorFra);
    v0AFra->AddNode(v0AFraR1,1);
    v0AFra->AddNode(v0AFraR2,1);
    v0AFra->AddNode(v0AFraR3,1);
    v0AFra->AddNode(v0AFraR4,1);
    v0AFra->AddNode(v0AFraR5,1);
    v0ASec->AddNode(v0AFra,1);

    /// Sensitive scintilator
    TGeoVolume *v0ASci = new TGeoVolumeAssembly("V0ASci");
    new TGeoTubeSeg( "sV0AR1b", fV0AR0+fV0AFraWd/2.,
		     fV0AR1-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR2b", fV0AR1+fV0AFraWd/2.,
		     fV0AR2-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR3b", fV0AR2+fV0AFraWd/2.,
		     fV0AR3-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR4b", fV0AR3+fV0AFraWd/2.,
		     fV0AR4-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    TGeoCompositeShape *sV0AR1 = new TGeoCompositeShape("sV0AR1","sV0AR1b-sV0ACha");
    TGeoCompositeShape *sV0AR2 = new TGeoCompositeShape("sV0AR2","sV0AR2b-sV0ACha");
    TGeoCompositeShape *sV0AR3 = new TGeoCompositeShape("sV0AR3","sV0AR3b-sV0ACha");
    TGeoCompositeShape *sV0AR4 = new TGeoCompositeShape("sV0AR4","sV0AR4b-sV0ACha");
    TGeoVolume *v0L1 = new TGeoVolume("V0L1",sV0AR1,medV0ASci);
    TGeoVolume *v0L2 = new TGeoVolume("V0L2",sV0AR2,medV0ASci);
    TGeoVolume *v0L3 = new TGeoVolume("V0L3",sV0AR3,medV0ASci);
    TGeoVolume *v0L4 = new TGeoVolume("V0L4",sV0AR4,medV0ASci);
    v0L1->SetLineColor(kV0AColorSci); v0L2->SetLineColor(kV0AColorSci);
    v0L3->SetLineColor(kV0AColorSci); v0L4->SetLineColor(kV0AColorSci);
    v0ASec->AddNode(v0L1,1);
    v0ASec->AddNode(v0L2,1);
    v0ASec->AddNode(v0L3,1);
    v0ASec->AddNode(v0L4,1);

    /// Non-sensitive scintilator
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR4;
      v0APts[1+8*i] = fV0AFraWd/2.;
      v0APts[2+8*i] = fV0AR4*sin45;
      v0APts[3+8*i] = (fV0AR4-fV0AFraWd)*sin45;
      v0APts[4+8*i] = fV0AR5/cos225*sin45+fV0AFraWd/2.*sin225;
      v0APts[5+8*i] = fV0AR5/cos225*sin45-fV0AFraWd/2.*cos225;
      v0APts[6+8*i] = fV0AR5/cos225-fV0AFraWd/2./ctg225;
      v0APts[7+8*i] = fV0AFraWd/2.;
    }
    new TGeoArb8("sV0AR5S1", fV0ASciWd/2., v0APts);
    new TGeoTubeSeg("sV0AR5S2", fV0AR4-(v0APts[6]-v0APts[0]),
		    fV0AR4+fV0AFraWd/2., fV0ASciWd/2.+2*preShape, 0, 45);
    TGeoCompositeShape *sV0AR5 = new TGeoCompositeShape("V0AR5","(sV0AR5S1 - sV0AR5S2)");
    TGeoVolume *v0AR5 = new TGeoVolume("V0AR5",sV0AR5,medV0ASci);
    v0AR5->SetLineColor(kV0AColorSci);
    v0ASci->AddNode(v0AR5,1);
    v0ASec->AddNode(v0ASci,1);

    /// Segment of innermost octagon
    TGeoVolume *v0ASup = new TGeoVolumeAssembly("V0ASup");
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = (fV0AR5-fV0AOctH1)/cos225;	v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = (fV0AR5-fV0AOctH1)/cos225*sin45;	v0APts[3+8*i] = (fV0AR5-fV0AOctH1)/cos225*sin45;
      v0APts[4+8*i] = fV0AR5/cos225*sin45;		v0APts[5+8*i] = fV0AR5/cos225*sin45;
      v0APts[6+8*i] = fV0AR5/cos225;			v0APts[7+8*i] = 0.;
    }
    TGeoArb8 *sV0AOct1 = new TGeoArb8("sV0AOct1", fV0AOctWd/2., v0APts);
    TGeoVolume *v0AOct1 = new TGeoVolume("V0AOct1",sV0AOct1,medV0ASup);
    v0AOct1->SetLineColor(kV0AColorOct);
    v0ASup->AddNode(v0AOct1,1,new TGeoTranslation(0,0,(fV0ASciWd+fV0AOctWd)/2.));
    v0ASup->AddNode(v0AOct1,2,new TGeoTranslation(0,0,-(fV0ASciWd+fV0AOctWd)/2.));

    /// Segment of outtermost octagon
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = (fV0AR6-fV0AOctH2)/cos225;	v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = (fV0AR6-fV0AOctH2)/cos225*sin45;  v0APts[3+8*i] = (fV0AR6-fV0AOctH2)/cos225*sin45;
      v0APts[4+8*i] = fV0AR6/cos225*sin45;		v0APts[5+8*i] = fV0AR6/cos225*sin45;
      v0APts[6+8*i] = fV0AR6/cos225;			v0APts[7+8*i] = 0.;
    }
    TGeoArb8 *sV0AOct2 = new TGeoArb8("sV0AOct2", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoVolume *v0AOct2 = new TGeoVolume("V0AOct2", sV0AOct2,medV0ASup);
    v0AOct2->SetLineColor(kV0AColorOct);
    v0ASup->AddNode(v0AOct2,1);
    v0ASec->AddNode(v0ASup,1);

    /// Bunch of fibers
    v0APts[ 0] = v0APts[ 2] = -12.5;
    v0APts[ 1] = v0APts[ 7] = (fV0ASciWd+fV0AOctWd)/2.-0.01;
    v0APts[ 3] = v0APts[ 5] = (fV0ASciWd+fV0AOctWd)/2.+0.01;
    v0APts[ 4] = v0APts[ 6] = +12.5;
    v0APts[ 8] = v0APts[10] = -0.5;
    v0APts[ 9] = v0APts[15] = 0.;
    v0APts[11] = v0APts[13] = 0.25;
    v0APts[12] = v0APts[14] = +0.5;
    TGeoArb8 *sV0AFib = new TGeoArb8("sV0AFib", (fV0AR6-fV0AR5-fV0AOctH2-0.006)/2., v0APts);
    TGeoVolume *v0AFib1 = new TGeoVolume("V0AFib1",sV0AFib,medV0AFib);
    TGeoVolume *v0AFib = new TGeoVolumeAssembly("V0AFib");
    TGeoRotation *rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateZ(-90.+22.5);
    v0AFib->AddNode(v0AFib1,1,rot);
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateY(180);
    rot->RotateZ(-90.+22.5);
    v0AFib->SetLineColor(kV0AColorFib);
    v0AFib->AddNode(v0AFib1,2,rot);
    v0ASec->AddNode(v0AFib,1,new TGeoTranslation((fV0AR6-fV0AOctH2+fV0AR5)*cos225/2.,
						 (fV0AR6-fV0AOctH2+fV0AR5)*sin225/2., 0));

    /// Plates
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0;			v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = fV0AR0*sin45;		v0APts[3+8*i] = fV0AR0*sin45;
      v0APts[4+8*i] = fV0AR6/cos225 * sin45;	v0APts[5+8*i] = fV0AR6/cos225*sin45;
      v0APts[6+8*i] = fV0AR6/cos225;		v0APts[7+8*i] = 0.;
    }
    TGeoArb8 *sV0APlaIn = new TGeoArb8("sV0APlaIn", (fV0APlaWd-2*fV0APlaAl)/2., v0APts);
    TGeoVolume *v0APlaIn = new TGeoVolume("V0APlaIn", sV0APlaIn, medV0APlaIn);
    TGeoArb8 *sV0APlaOu = new TGeoArb8("sV0APlaOu", fV0APlaAl/2., v0APts);
    TGeoVolume *v0APlaOu = new TGeoVolume("V0APlaOu", sV0APlaOu, medV0APlaOu);
    v0APlaIn->SetLineColor(kV0AColorPlaIn); v0APlaOu->SetLineColor(kV0AColorPlaOu);
    TGeoVolume *v0APla = new TGeoVolumeAssembly("V0APla");
    v0APla->AddNode(v0APlaIn,1);
    v0APla->AddNode(v0APlaOu,1,new TGeoTranslation(0,0,(fV0APlaWd-fV0APlaAl)/2.));
    v0APla->AddNode(v0APlaOu,2,new TGeoTranslation(0,0,-(fV0APlaWd-fV0APlaAl)/2.));
    v0ASec->AddNode(v0APla,1,new TGeoTranslation(0,0,(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    v0ASec->AddNode(v0APla,2,new TGeoTranslation(0,0,-(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));

    /// PMBox
    TGeoVolume* v0APM = new TGeoVolumeAssembly("V0APM");
    new TGeoBBox("sV0APMB1", fV0APMBWd/2., fV0APMBHt/2., fV0APMBTh/2.);
    new TGeoBBox("sV0APMB2", fV0APMBWd/2.-fV0APMBWdW, fV0APMBHt/2.-fV0APMBHtW, fV0APMBTh/2.-fV0APMBThW);
    TGeoCompositeShape *sV0APMB = new TGeoCompositeShape("sV0APMB","sV0APMB1-sV0APMB2");
    TGeoVolume *v0APMB = new TGeoVolume("V0APMB",sV0APMB, medV0APMAlum);
    v0APMB->SetLineColor(kV0AColorPMA);
    v0APM->AddNode(v0APMB,1);

    /// PMTubes
    TGeoTube *sV0APMT1 = new TGeoTube("sV0APMT1", fV0APMTR1, fV0APMTR2, fV0APMTH/2.);
    TGeoVolume *v0APMT1 = new TGeoVolume("V0APMT1", sV0APMT1, medV0APMGlass);
    TGeoTube *sV0APMT2 = new TGeoTube("sV0APMT2", fV0APMTR3, fV0APMTR4, fV0APMTH/2.);
    TGeoVolume *v0APMT2 = new TGeoVolume("V0APMT2", sV0APMT2, medV0APMAlum);
    TGeoVolume *v0APMT = new TGeoVolumeAssembly("V0APMT");
    TGeoTube *sV0APMTT = new TGeoTube("sV0APMTT", 0., fV0APMTR4, fV0APMTB/2.);
    TGeoVolume *v0APMTT = new TGeoVolume("V0APMT1", sV0APMTT, medV0APMAlum);
    v0APMT1->SetLineColor(kV0AColorPMG);
    v0APMT2->SetLineColor(kV0AColorPMA);
    v0APMTT->SetLineColor(kV0AColorPMA);
    rot = new TGeoRotation("rot", 90, 0, 180, 0, 90, 90);
    v0APMT->AddNode(v0APMT1,1,rot);
    v0APMT->AddNode(v0APMT2,1,rot);
    v0APMT->AddNode(v0APMTT,1,new TGeoCombiTrans(0,-(fV0APMTH+fV0APMTB)/2.,0,rot));
    double autoShift = (fV0APMBWd-2*fV0APMBWdW)/4.;
    v0APM->AddNode(v0APMT, 1, new TGeoTranslation(-1.5*autoShift, 0, 0));
    v0APM->AddNode(v0APMT, 2, new TGeoTranslation(-0.5*autoShift, 0, 0));
    v0APM->AddNode(v0APMT, 3, new TGeoTranslation(+0.5*autoShift, 0, 0));
    v0APM->AddNode(v0APMT, 4, new TGeoTranslation(+1.5*autoShift, 0, 0));

    /// PM
    rot = new TGeoRotation("rot");
    rot->RotateX(90-fV0APMBAng);
    rot->RotateZ(-90.+22.5);
    double cosAngPMB = TMath::Cos(fV0APMBAng*TMath::DegToRad());
    double sinAngPMB = TMath::Sin(fV0APMBAng*TMath::DegToRad());
    double shiftZ = fV0APMBHt/2. * cosAngPMB
      -   ( fV0ASciWd + 2 * fV0AOctWd + 2 * fV0APlaWd )/2.   -   fV0APMBTh/2. * sinAngPMB;
    double shiftR = fV0AR6  +  fV0APMBHt/2. * sinAngPMB  +  fV0APMBTh/2. * cosAngPMB;
    v0ASec->AddNode(v0APM,1, new TGeoCombiTrans( shiftR*cos225, shiftR*sin225, shiftZ, rot));

    /// End of sector definition
    ////////////////////////////

    /// Replicate sectors
    TGeoVolume *v0LE = new TGeoVolumeAssembly("V0LE");
    for(int i=0; i<8; i++) {
      TGeoRotation *rot = new TGeoRotation("rot", 90., i*45.+90, 90., 90.+i*45.+90, 0., 0.);
      v0LE->AddNode(v0ASec,i+1,rot);  /// modificacion +1 anhadido
    }
  
    /// Basis Construction
    rot = new TGeoRotation("rot"); rot->RotateX(90-fV0APMBAng); rot->RotateZ(-22.5);
    TGeoCombiTrans *pos1 = new TGeoCombiTrans("pos1", shiftR*sin225, shiftR*cos225, shiftZ, rot);
    pos1->RegisterYourself();
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR6/cos225*sin45;  v0APts[1+8*i] = fV0AR6/cos225*sin45;
      v0APts[2+8*i] = 0;		    v0APts[3+8*i] = fV0AR6/cos225;
      v0APts[4+8*i] = 0;		    v0APts[5+8*i] = fV0AR6/cos225+fV0APlaEx;
      v0APts[6+8*i] = fV0AR6/cos225-(fV0AR6/cos225+fV0APlaEx)/ctg225;
      v0APts[7+8*i] = fV0AR6/cos225+fV0APlaEx;
    }
    new TGeoArb8("sV0APlaExIn1", (fV0APlaWd-2*fV0APlaAl)/2., v0APts);
    new TGeoArb8("sV0APlaExOu1", fV0APlaAl/2., v0APts);
    TGeoCompositeShape *sV0APlaExIn = new TGeoCompositeShape("sV0APlaExIn","sV0APlaExIn1-sV0APMB1:pos1");
    TGeoVolume *v0APlaExIn = new TGeoVolume("V0APlaExIn", sV0APlaExIn, medV0APlaIn);
    TGeoCompositeShape *sV0APlaExOu = new TGeoCompositeShape("sV0APlaExOu","sV0APlaExOu1-sV0APMB1:pos1");
    TGeoVolume *v0APlaExOu = new TGeoVolume("V0APlaExOu", sV0APlaExOu, medV0APlaOu);
    v0APlaExIn->SetLineColor(kV0AColorPlaIn); v0APlaExOu->SetLineColor(kV0AColorPlaOu);
    TGeoVolume *v0APlaEx = new TGeoVolumeAssembly("V0APlaEx");
    v0APlaEx->AddNode(v0APlaExIn,1);
    v0APlaEx->AddNode(v0APlaExOu,1,new TGeoTranslation(0,0,(fV0APlaWd-fV0APlaAl)/2.));
    v0APlaEx->AddNode(v0APlaExOu,2,new TGeoTranslation(0,0,-(fV0APlaWd-fV0APlaAl)/2.));
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR6/cos225-(fV0AR6/cos225+fV0APlaEx)/ctg225-fV0ABasHt*sin45;
      v0APts[1+8*i] = fV0AR6/cos225+fV0APlaEx-fV0ABasHt*sin45;
      v0APts[2+8*i] = 0;  v0APts[3+8*i] = fV0AR6/cos225+fV0APlaEx-fV0ABasHt;
      v0APts[4+8*i] = 0;  v0APts[5+8*i] = fV0AR6/cos225+fV0APlaEx;
      v0APts[6+8*i] = fV0AR6/cos225-(fV0AR6/cos225+fV0APlaEx)/ctg225;
      v0APts[7+8*i] = fV0AR6/cos225+fV0APlaEx;
    }
    new TGeoArb8("sV0ABas1", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoCompositeShape *sV0ABas = new TGeoCompositeShape("sV0ABas","sV0ABas1-sV0APMB1:pos1");
    TGeoVolume *v0ABas = new TGeoVolume("V0ABas", sV0ABas, medV0ABas);
    v0ABas->SetLineColor(kV0AColorBas);
    TGeoVolume *v0ABasis = new TGeoVolumeAssembly("V0ABasis");
    rot = new TGeoRotation("rot",90.,180.,90.,90.,0.,0.);
    v0ABasis->AddNode(v0APlaEx,1, new TGeoTranslation(0,0,(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    v0ABasis->AddNode(v0APlaEx,2, new TGeoTranslation(0,0,-(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    v0ABasis->AddNode(v0APlaEx,3, new TGeoCombiTrans(0,0,(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.,rot));
    v0ABasis->AddNode(v0APlaEx,4, new TGeoCombiTrans(0,0,-(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.,rot));
    v0ABasis->AddNode(v0ABas,1);
    v0ABasis->AddNode(v0ABas,2,rot);
    rot = new TGeoRotation("rot");
    rot->RotateZ(180);
    v0LE->AddNode(v0ABasis,1,rot);

  // Adding detectors to top volume
    TGeoVolume *vZERO = new TGeoVolumeAssembly("VZERO");
    vZERO->AddNode(v0RI,1,new TGeoTranslation(0, 0, -zdet));
    vZERO->AddNode(v0LE,1,new TGeoTranslation(0, 0, +340));
    top->AddNode(vZERO,1);
}

//_____________________________________________________________________________
void AliVZEROv7::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  // 
  TString vpC = "/ALIC_1/VZERO_1/V0RI_1";
  TString vpA = "/ALIC_1/VZERO_1/V0LE_1";
  TString snC = "VZERO/V0C";
  TString snA = "VZERO/V0A";
  
  if(!gGeoManager->SetAlignableEntry(snC.Data(),vpC.Data()))
    AliFatal(Form("Alignable entry %s not created. Volume path %s not valid", snC.Data(),vpC.Data()));
  if(!gGeoManager->SetAlignableEntry(snA.Data(),vpA.Data()))
    AliFatal(Form("Alignable entry %s not created. Volume path %s not valid", snA.Data(),vpA.Data()));

} 

//_____________________________________________________________________________
void AliVZEROv7::CreateMaterials()
{

// Creates materials used for geometry 

  AliDebug(2,"Create materials");
  // Parameters for simulation scope
  Int_t     fieldType       = gAlice->Field()->Integ();     // Field type 
  Double_t  maxField        = gAlice->Field()->Max();       // Field max.
  Double_t  maxBending      = 10;    // Max Angle
  Double_t  maxStepSize     = 0.01;  // Max step size 
  Double_t  maxEnergyLoss   = 1;     // Max Delta E
  Double_t  precision       = 0.003; // Precision
  Double_t  minStepSize     = 0.003; // Minimum step size 

  Int_t    id;
  Double_t a, z, radLength, absLength;
  Float_t density, as[4], zs[4], ws[4];

// Parameters  for V0CPrePlates: Aluminium
   a = 26.98; 
   z = 13.00;
   density     = 2.7;
   radLength   = 8.9;
   absLength   = 37.2;
   id = 2;
   AliMaterial( id, "V0CAlu", a, z, density, radLength, absLength, 0, 0);
   AliMedium(id, "V0CAlu", id, 1, fieldType, maxField, maxBending, maxStepSize,
	     maxEnergyLoss, precision, minStepSize);
		    
// Parameters  for V0CPlates: Carbon 
   a = 12.01; 
   z =  6.00;
   density   = 2.265;
   radLength = 18.8;
   absLength = 49.9;
   id = 3;
   AliMaterial(id, "V0CCar",  a, z, density, radLength, absLength, 0, 0);
   AliMedium(id, "V0CCar", id, 1, fieldType, maxField, maxBending, maxStepSize,
	     maxEnergyLoss, precision, minStepSize);
	    
// Parameters  for V0Cscintillator: BC408
   as[0] = 1.00794;	as[1] = 12.011;
   zs[0] = 1.;		zs[1] = 6.;
   ws[0] = 1.;	        ws[1] = 1.;
   density      = 1.032;
   id           = 4;
   AliMixture(id, "V0CSci", as, zs, density, -2, ws);
   AliMedium(id,"V0CSci", id, 1, fieldType, maxField, maxBending, maxStepSize,
	     maxEnergyLoss, precision, minStepSize);

// Parameters for V0Ascintilator: BC404
   as[0] = 1.00794;	as[1] = 12.011;
   zs[0] = 1.;		zs[1] = 6.;
   ws[0] = 5.21;	ws[1] = 4.74;
   density      = 1.032;
   id           = 5;
   AliMixture(id, "V0ASci", as, zs, density, -2, ws);
   AliMedium(id,  "V0ASci", id, 1, fieldType, maxField, maxBending, maxStepSize,
	     maxEnergyLoss, precision, minStepSize);

// Parameters for V0ALuc: Lucita but for the simulation BC404
   as[0] = 1.00794;	as[1] = 12.011;
   zs[0] = 1.;		zs[1] = 6.;
   ws[0] = 5.21;	ws[1] = 4.74;
   density      = 1.032;
   id           = 6;
   AliMixture(id, "V0ALuc", as, zs, density, -2, ws);
   AliMedium(id, "V0ALuc", id, 1, fieldType, maxField, maxBending, maxStepSize,
	     maxEnergyLoss, precision, minStepSize);

// Parameters for V0Aplate: EuroComposite - EC-PI 626 PS - AlMg3
   as[0] = 26.982;	as[1] = 24.305;
   zs[0] = 13.;		zs[1] = 12.;
   ws[0] = 1.;		ws[1] = 3.;
   density      = 3.034;
   id           = 7;
   AliMixture(id, "V0APlaOu", as, zs, density, -2, ws);
   AliMedium(id, "V0APlaOu", id, 1, fieldType, maxField, maxBending, maxStepSize,
	     maxEnergyLoss, precision, minStepSize);

// Parameters for V0Aplate: EuroComposite - EC-PI 626 PS - EC-PI 6.4-42
   as[0] = 1.00794;	as[1] = 12.011;
   zs[0] = 1.;		zs[1] = 6.;
   ws[0] = 5.21;	ws[1] = 4.74;
   density      = 0.042;
   id           = 8;
   AliMixture(id, "V0APlaIn", as, zs, density, -2, ws);
   AliMedium(id, "V0APlaIn", id, 1, fieldType, maxField, maxBending, maxStepSize,
	     maxEnergyLoss, precision, minStepSize);

// Parameters for V0Afiber: BC9929AMC Plastic Scintillating Fiber from Saint-Gobain
   as[0] = 1.00794;	as[1] = 12.011;
   zs[0] = 1.;		zs[1] = 6.;
   ws[0] = 4.82;	ws[1] = 4.85;
   density      = 1.05;
   id           = 9;
   AliMixture(id, "V0AFib", as, zs, density, -2, ws);
   AliMedium(id, "V0AFib", id, 1, fieldType, maxField, maxBending, maxStepSize,
	     maxEnergyLoss, precision, minStepSize);

// Parameters for V0APMA: Aluminium
   a = 26.98; 
   z = 13.00;
   density     = 2.7;
   radLength   = 8.9;
   absLength   = 37.2;
   id = 10;
   AliMaterial(id, "V0APMA",  a, z, density, radLength, absLength, 0, 0);
   AliMedium(id, "V0APMA", id, 1, fieldType, maxField, maxBending, maxStepSize,
	     maxEnergyLoss, precision, minStepSize);

// Parameters for V0APMG: Glass for the simulation Aluminium
   a = 26.98; 
   z = 13.00;
   density   = 2.7;
   radLength = 8.9;
   absLength = 37.2;
   id = 11;
   AliMaterial(id, "V0APMG",  a, z, density, radLength, absLength, 0, 0);
   AliMedium(id, "V0APMG", id, 1, fieldType, maxField, maxBending, maxStepSize,
	     maxEnergyLoss, precision, minStepSize);
}

//_____________________________________________________________________________
void AliVZEROv7::DrawModule() const
{
//  Drawing is done in DrawVZERO.C

   AliDebug(2,"DrawModule");
}


//_____________________________________________________________________________
void AliVZEROv7::DrawGeometry() 
{
//  Drawing of V0 geometry done in DrawV0.C

   AliDebug(2,"DrawGeometry");
}

//_____________________________________________________________________________
void AliVZEROv7::Init()
{
// Initialises version of the VZERO Detector given in Config
// Just prints an information message

//   AliInfo(Form("VZERO version %d initialized \n",IsVersion()));
   
   AliDebug(1,"VZERO version 7 initialized");
   AliVZERO::Init();  
}

//_____________________________________________________________________________
void AliVZEROv7::StepManager()
{
  // Step Manager, called at each step 

  Int_t     copy;
  static    Int_t   vol[4];
  static    Float_t hits[21];
  static    Float_t eloss, tlength;
  static    Int_t   nPhotonsInStep;
  static    Int_t   nPhotons; 
  static    Int_t   numStep;
  Int_t     ringNumber;
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
  bool   hitOnV0C = true;
  double lightYield;
  double lightAttenuation;
  double nMeters;
  double fibToPhot;
  if      ( gMC->CurrentVolID(copy) == idV0R1 || gMC->CurrentVolID(copy) == idV0L1 )
    ringNumber = 1;
  else if ( gMC->CurrentVolID(copy) == idV0R2 || gMC->CurrentVolID(copy) == idV0L2 )
    ringNumber = 2;  
  else if ( gMC->CurrentVolID(copy) == idV0R3 || gMC->CurrentVolID(copy) == idV0R4
	    || gMC->CurrentVolID(copy) == idV0L3 ) ringNumber = 3;
  else if ( gMC->CurrentVolID(copy) == idV0R5 || gMC->CurrentVolID(copy) == idV0R6
	    || gMC->CurrentVolID(copy) == idV0L4 ) ringNumber = 4;	       
  else ringNumber = 0;
  if  (ringNumber) {
    if (gMC->CurrentVolID(copy) == idV0L1 || gMC->CurrentVolID(copy) == idV0L2 ||
	gMC->CurrentVolID(copy) == idV0L3 || gMC->CurrentVolID(copy) == idV0L4)
      hitOnV0C = false;
    destep = gMC->Edep();
    step   = gMC->TrackStep();
    if (hitOnV0C) {
      lightYield = fV0CLightYield;
      lightAttenuation = fV0CLightAttenuation;
      nMeters = fV0CnMeters;
      fibToPhot = fV0CFibToPhot;
    } else {
      lightYield = fV0ALightYield;
      lightAttenuation = fV0ALightAttenuation;
      nMeters = fV0AnMeters;
      fibToPhot = fV0AFibToPhot;
    }
    nPhotonsInStep  = Int_t(destep / (lightYield *1e-9) );	
    nPhotonsInStep  = gRandom->Poisson(nPhotonsInStep);
    eloss    += destep;
    tlength  += step; 	 
    if ( gMC->IsTrackEntering() ) { 
      nPhotons  =  nPhotonsInStep;
      gMC->TrackPosition(fTrackPosition);
      gMC->TrackMomentum(fTrackMomentum);
      Float_t pt  = TMath::Sqrt( fTrackMomentum.Px() * fTrackMomentum.Px()
				 + fTrackMomentum.Py() * fTrackMomentum.Py() );
      TParticle *par = gAlice->GetMCApp()->Particle(gAlice->GetMCApp()->GetCurrentTrackNumber());
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
      hits[14] = par->Vx();
      hits[15] = par->Vy();
      hits[16] = par->Vz();
      tlength  = 0.0;
      eloss    = 0.0;	    

      //////////////////////////
      ///// Display V0A geometry
      //      if (!hitOnV0C) {
      //      	FILE *of;
      //      	of = fopen("V0A.out", "a");
      //      	// x, y, z, ringnumber, cellid
      //      	fprintf( of, "%f %f %f %f %d \n",  hits[0], hits[1], hits[2], hits[8], GetCellId (vol, hits) );
      //      	fclose(of);
      //      }
      //////////////////////////
    }
    nPhotons  = nPhotons + nPhotonsInStep;
    if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()){
      nPhotons = nPhotons - Int_t((Float_t(nPhotons) * lightAttenuation * nMeters));
      nPhotons = nPhotons - Int_t( Float_t(nPhotons) * fibToPhot);
      hits[17] = eloss;
      hits[18] = tlength;
      hits[19] = nPhotons;
      hits[20] = GetCellId (vol, hits);
      AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);
      tlength         = 0.0;
      eloss           = 0.0; 
      nPhotons        = 0;
      nPhotonsInStep  = 0;
      numStep         = 0;  
    }
  }
}

//_____________________________________________________________________________
void AliVZEROv7::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
//  Adds a VZERO hit

  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliVZEROhit(fIshunt,track,vol,hits);
}

//_____________________________________________________________________________
void AliVZEROv7::AddDigits(Int_t *tracks, Int_t* digits) 
{
//  Adds a VZERO digit

   TClonesArray  &ldigits = *fDigits;
   new(ldigits[fNdigits++]) AliVZEROdigit(tracks, digits);
}

//_____________________________________________________________________________
void AliVZEROv7::MakeBranch(Option_t *option)
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
Int_t AliVZEROv7::GetCellId(Int_t *vol, Float_t *hits) 
{
  //   Returns Id of scintillator cell
  //   Right side from  0 to 47 
  //   Left  side from 48 to 79
  //   hits[8] = ring number (1 to 4)
  //   vol[1]  = copy number (1 to 8)

  Int_t index      = vol[1];
  Int_t ringNumber = Int_t(hits[8]);
  fCellId          = 0;

  Float_t phi = Float_t(TMath::ATan2(Double_t(hits[1]),Double_t(hits[0])) ); 
  Float_t kRaddeg = 180.0/TMath::Pi();
  phi = kRaddeg * phi;

  if (index < 7) index = index + 8;

  if (hits[2] < 0.0) {
    if(ringNumber < 3) {
      index = (index - 7) + ( ( ringNumber - 1 ) * 8);
    } else if (ringNumber >= 3) { 
      if ( gMC->CurrentVolID(vol[1]) == gMC->VolId("V0R3") || gMC->CurrentVolID(vol[1])
	   == gMC->VolId("V0R5") )  index = (index*2-14)+((ringNumber-2)*16);
      if ( gMC->CurrentVolID(vol[1]) == gMC->VolId("V0R4") || gMC->CurrentVolID(vol[1])
	   == gMC->VolId("V0R6") )  index = (index*2-13)+((ringNumber-2)*16);
    }
    fCellId   = index;           
  } else if (hits[2] > 0.0) {
    index = (index - 7 + 48) + ( ( ringNumber - 1 ) * 8);
    fCellId   = index;
  }

  return fCellId;
}
