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
//   V0L (now V0A) sits between Z values +325.0 and +330.0 cm        // 
//   New coordinate system has been implemented in october 2003      //
//   Revision of the V0A part by Lizardo Valencia  in July 2008      //
//                                                                   //
/////////////////////////////////////////////////////////////////////// 

// --- Standard libraries ---
#include <Riostream.h>

// --- ROOT libraries ---
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include <TMath.h>
#include <TParticle.h>
#include <TVirtualMC.h>

#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include "TGeoTube.h"
#include "TGeoArb8.h"
#include "TGeoCompositeShape.h"
#include <TTree.h>

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliVZEROLoader.h"
#include "AliVZEROdigit.h"
#include "AliVZEROhit.h"
#include "AliVZEROv7.h"
#include "AliLog.h"
#include "AliTrackReference.h"
 
ClassImp(AliVZEROv7)

//_____________________________________________________________________________
AliVZEROv7:: AliVZEROv7():AliVZERO(),
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
   fV0COffsetFibers(1.125),
   fV0CLightYield(93.75),
   fV0CLightAttenuation(0.05),
   fV0CnMeters(15.0),
   fV0CFibToPhot(0.3),
   fV0AR0(4.2),
   fV0AR1(7.6), 
   fV0AR2(13.8), 
   fV0AR3(22.7),
   fV0AR4(41.3), 
   fV0AR5(43.3), 
   fV0AR6(72.6),
   fV0AR7(92.0), // Distance from origin to outtermost intersection sector7 and sector8 
   fV0ASciWd(2.5), 
   fV0APlaWd(0.5), 
   fV0APlaAl(0.06), 
   fV0AOctWd(0.75), 
   fV0AFraWd(0.2),
   fV0AOctH1(1.0), 
   fV0AOctH2(2.0), 
   fV0ABasHt(2.0),
   fV0AFibRd(0.1),
   fV0APlaEx(4.4),
   fV0APMBWd(24.6), 
   fV0APMBHt(22.0), 
   fV0APMBTh(7.1), 
   fV0APMBWdW(0.3), 
   fV0APMBHtW(1.0),
   fV0APMBAng(30.0), 
   fV0APMBThW(0.3), 
   fV0APMTR1(2.44), 
   fV0APMTR2(2.54), 
   fV0APMTR3(2.54),
   fV0APMTR4(2.70), 
   fV0APMTH(10.0), 
   fV0APMTB(1.0),
   fV0AFEEBWd(26.5),
   fV0AFEEBHt(20.5),
   fV0AFEEBTh(7.5),
   fV0AnMeters(fV0AR6*0.01),
   fV0ALightYield(93.75),
   fV0ALightAttenuation(0.05), 
   fV0AFibToPhot(0.3),
   fVersion(7)
{
// Standard default constructor 
}

//_____________________________________________________________________________
AliVZEROv7::AliVZEROv7(const char *name, const char *title):AliVZERO(name,title),
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
   fV0COffsetFibers(1.125),
   fV0CLightYield(93.75),
   fV0CLightAttenuation(0.05),
   fV0CnMeters(15.0),
   fV0CFibToPhot(0.3),
   fV0AR0(4.2),
   fV0AR1(7.6), 
   fV0AR2(13.8), 
   fV0AR3(22.7),
   fV0AR4(41.3), 
   fV0AR5(43.3), 
   fV0AR6(72.6),
   fV0AR7(92.0), // Distance from origin to outtermost intersection of sector7 and sector8 
   fV0ASciWd(2.5), 
   fV0APlaWd(0.5), 
   fV0APlaAl(0.06), 
   fV0AOctWd(0.75), 
   fV0AFraWd(0.2),
   fV0AOctH1(1.0), 
   fV0AOctH2(2.0), 
   fV0ABasHt(2.0),
   fV0AFibRd(0.1),
   fV0APlaEx(4.4),
   fV0APMBWd(24.6), 
   fV0APMBHt(22.0), 
   fV0APMBTh(7.1), 
   fV0APMBWdW(0.3), 
   fV0APMBHtW(1.0),
   fV0APMBAng(30.0), 
   fV0APMBThW(0.3), 
   fV0APMTR1(2.44), 
   fV0APMTR2(2.54), 
   fV0APMTR3(2.54),
   fV0APMTR4(2.70), 
   fV0APMTH(10.0), 
   fV0APMTB(1.0),
   fV0AFEEBWd(26.5),
   fV0AFEEBHt(20.5),
   fV0AFEEBTh(7.5),		   
   fV0AnMeters(fV0AR6*0.01),
   fV0ALightYield(93.75),
   fV0ALightAttenuation(0.05),
   fV0AFibToPhot(0.3),
   fVersion(7)


{
// Standard constructor for V-zero Detector  version 7

  AliDebug(2,"Create VZERO object ");

//  fVersion            =     7;  // version number

//   // V0C Parameters related to geometry: All in cm
//   fV0CHeight1         =    2.5; // height of cell 1
//   fV0CHeight2         =    4.4; // height of cell 2
//   fV0CHeight3         =    7.4; // height of cell 3
//   fV0CHeight4         =   12.5; // height of cell 4
//   fV0CRMin            =    4.6; // inner radius of box
//   fV0CRBox            =   38.0; // outer radius of box
//   fV0CLidThickness    =   0.30; // thickness of Carbon lid
//   fV0CCellThickness   =   2.00; // thickness of elementary cell
//   fV0CBoxThickness    =   4.70; // thickness of V0C Box
//   fV0COffsetFibers    =    1.0; // offset to output fibers
//   // V0C Parameters related to light output
//   fV0CLightYield         =  93.75; // Light yield in BC408 (93.75 eV per photon)
//   fV0CLightAttenuation   =   0.05; // Light attenuation in fiber (0.05 per meter)
//   fV0CnMeters            =   15.0; // Number of meters of clear fibers to PM
//   fV0CFibToPhot          =    0.3; // Attenuation at fiber-photocathode interface
// 
//   // V0A Parameters related to geometry: All in cm
//   fV0AR0     =  4.2;  // Radius of hole
//   fV0AR1     =  7.6;  // Maximun radius of 1st cell
//   fV0AR2     = 13.8; // Maximun radius of 2nd cell
//   fV0AR3     = 22.7; // Maximun radius of 3rd cell
//   fV0AR4     = 41.3; // Maximun radius of 4th cell
//   fV0AR5     = 43.3; // Radius circunscrite to innermost octagon
//   fV0AR6     = 68.0; // Radius circunscrite to outtermost octagon
//   fV0ASciWd  =  2.5;  // Scintillator thickness 
//   fV0APlaWd  =  0.5;  // Plates thinckness
//   fV0APlaAl  = 0.06; // Plates AlMg3 thinckness
//   fV0AOctWd  = 0.75; // Innermost octagon thickness
//   fV0AOctH1  =  1.0;  // Height of innermost octagon
//   fV0AOctH2  =  2.0;  // Height of outtermost octagon
//   fV0AFibRd  =  0.1;  // Radius of Fiber
//   fV0AFraWd  =  0.2;  // Support Frame thickness
//   fV0APMBWd  = 24.6;  // Width of PM Box
//   fV0APMBHt  = 22.0;  // Height of PM Box
//   fV0APMBTh  =  7.1;  // Thickness of PM Box
//   fV0APMBWdW =  0.3;  // Thickness of PM Box Side1 Wall
//   fV0APMBHtW =  1.0;  // Thickness of PM Box Side2 Wall
//   fV0APMBThW =  0.3;  // Thickness of PM Box Top Wall
//   fV0APMBAng = 30.0;  // Angle between PM Box and Support
//   fV0APMTR1  = 2.44;  // PMT Glass
//   fV0APMTR2  = 2.54;  // PMT Glass
//   fV0APMTR3  = 2.54;  // PMT Cover
//   fV0APMTR4  = 2.70;  // PMT Cover
//   fV0APMTH   = 10.0;  // PMT Height
//   fV0APMTB   =  1.0;  // PMT Basis
//   fV0APlaEx  =  4.4;  // Plates Extension height
//   fV0ABasHt  =  2.0;  // Basis Height
//   // V0A Parameters related to light output
//   fV0ALightYield         =  93.75;      // Light yield in BC404
//   fV0ALightAttenuation   =   0.05;      // Light attenuation in WLS fiber, per meter
//   fV0AnMeters            = fV0AR6*0.01; // Tentative value, in meters
//   fV0AFibToPhot          =    0.3;      // Attenuation at fiber-photocathode interface
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
    //    TGeoTube   *sV0CA = new TGeoTube("V0CA", partube[0], partube[1], partube[2]);
    //    TGeoVolume *v0CA  = new TGeoVolume("V0CA",sV0CA,medV0CCar);
    //    TGeoTranslation *tr2 = new TGeoTranslation(0.,0., fV0CBoxThickness/2.0-partube[2]);
    //    TGeoTranslation *tr3 = new TGeoTranslation(0.,0.,-fV0CBoxThickness/2.0+partube[2]);
    //    v0RI->AddNode(v0CA,1,tr2);
    //    v0RI->AddNode(v0CA,2,tr3);
    //    v0CA->SetLineColor(kYellow);

    Float_t rInt1 = 11.5, rOut1 = 20.0, rInt2 = 9.0;

    TGeoTube   *sV0CA4 = new TGeoTube("V0CA4", partube[0], rInt2, partube[2] - 0.1);
    TGeoVolume *v0CA4  = new TGeoVolume("V0CA4",sV0CA4,medV0CCar);
    TGeoTranslation *tr21 = new TGeoTranslation(0.,0., fV0CBoxThickness/2.0-partube[2] + 0.1);
    v0RI->AddNode(v0CA4,1,tr21);
    v0CA4->SetLineColor(kYellow);

    TGeoTube   *sV0CA5 = new TGeoTube("V0CA5", rInt2, partube[1], partube[2]);
    TGeoVolume *v0CA5  = new TGeoVolume("V0CA5",sV0CA5,medV0CCar);
    TGeoTranslation *tr22 = new TGeoTranslation(0.,0., fV0CBoxThickness/2.0-partube[2]);
    v0RI->AddNode(v0CA5,1,tr22);
    v0CA5->SetLineColor(kYellow);

    TGeoTube   *sV0CA1 = new TGeoTube("V0CA1",partube[0], rInt1, partube[2]);
    TGeoVolume *v0CA1  = new TGeoVolume("V0CA1",sV0CA1,medV0CCar);
    TGeoTranslation *tr31 = new TGeoTranslation(0.,0.,-fV0CBoxThickness/2.0+partube[2]);
    v0RI->AddNode(v0CA1,1,tr31);
    v0CA1->SetLineColor(kYellow);

    TGeoTube   *sV0CA2 = new TGeoTube("V0CA2", rInt1, rOut1, partube[2] - 0.1);
    TGeoVolume *v0CA2  = new TGeoVolume("V0CA2",sV0CA2,medV0CCar);
    TGeoTranslation *tr32 = new TGeoTranslation(0.,0.,-fV0CBoxThickness/2.0+partube[2] - 0.1);
    v0RI->AddNode(v0CA2,1,tr32);
    v0CA2->SetLineColor(kYellow);

    TGeoTube   *sV0CA3 = new TGeoTube("V0CA3", rOut1, partube[1], partube[2]);
    TGeoVolume *v0CA3  = new TGeoVolume("V0CA3",sV0CA3,medV0CCar);
    TGeoTranslation *tr33 = new TGeoTranslation(0.,0.,-fV0CBoxThickness/2.0+partube[2]);
    v0RI->AddNode(v0CA3,1,tr33);
    v0CA3->SetLineColor(kYellow);

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
    TGeoTranslation *tr4 = new TGeoTranslation(0.,0., offset);
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
    TGeoTranslation *tr5 = new TGeoTranslation(0.0,0.2, offset - fV0COffsetFibers);
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
    TGeoTranslation *tr6 = new TGeoTranslation(0.,0.2, offset - 2.0*fV0COffsetFibers);
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
    TGeoTranslation *tr7 = new TGeoTranslation(0.,0.0, offset - 2.0*fV0COffsetFibers + 0.25);					      
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
  // Revision by Lizardo VALENCIA, UNAM Mexico in July 2008

    const int kV0AColorSci   = 5;
    const int kV0AColorPlaIn = 3;
    const int kV0AColorPlaOu = 41;
    const int kV0AColorOct   = 7;
    const int kV0AColorFra   = 6;
    const int kV0AColorFib   = 11;
    const int kV0AColorPMG   = 1;
    const int kV0AColorPMA   = 2;
    const int kV0AColorFibGlass = 4; 
    TGeoMedium *medV0ASci = gGeoManager->GetMedium("VZERO_V0ASci");
    TGeoMedium *medV0APlaIn = gGeoManager->GetMedium("VZERO_V0APlaIn");
    TGeoMedium *medV0APlaOu = gGeoManager->GetMedium("VZERO_V0APlaOu");
    TGeoMedium *medV0ASup = gGeoManager->GetMedium("VZERO_V0APMA");
    TGeoMedium *medV0AFra = gGeoManager->GetMedium("VZERO_V0ALuc");
    TGeoMedium *medV0AFib = gGeoManager->GetMedium("VZERO_V0AFib");
    TGeoMedium *medV0APMGlass = gGeoManager->GetMedium("VZERO_V0APMG");
    TGeoMedium *medV0APMAlum = gGeoManager->GetMedium("VZERO_V0APMA");
    TGeoMedium *medV0AFibGlass = gGeoManager->GetMedium("VZERO_V0AFibGlass");
    double pi = TMath::Pi();
    double sin225   = TMath::Sin(pi/8.);
    double cos225   = TMath::Cos(pi/8.);
    double sin45    = TMath::Sin(pi/4.); // lucky: Sin45=Cos45
    double cos45    = TMath::Cos(pi/4.); 
    double v0APts[16];
    double sin654   = TMath::Sin(1.14);
    double cos654   = TMath::Cos(1.14);
    
    //Defining the master volume for V0A
    TGeoVolume *v0LE = new TGeoVolumeAssembly("V0LE");

    /// Definition sector 1
    TGeoVolume *v0ASec1 = new TGeoVolumeAssembly("V0ASec1");
        
    /// For boolean sustraction
    double preShapeSec1 = 0.2;
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0-fV0AFraWd/2.-preShapeSec1;  v0APts[1+8*i] = -preShapeSec1;
      v0APts[2+8*i] = fV0AR0-fV0AFraWd/2.-preShapeSec1;  v0APts[3+8*i] = fV0AFraWd/2.;
      v0APts[4+8*i] = fV0AR4+fV0AFraWd/2.+preShapeSec1;  v0APts[5+8*i] = fV0AFraWd/2.;
      v0APts[6+8*i] = fV0AR4+fV0AFraWd/2.+preShapeSec1;  v0APts[7+8*i] = -preShapeSec1;
    }
    new TGeoArb8("sV0ACha1Sec1",fV0ASciWd/1.5,v0APts);
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0*sin45-preShapeSec1;
      v0APts[1+8*i] = (fV0AR0-fV0AFraWd)*sin45-preShapeSec1;
      v0APts[2+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45-preShapeSec1;
      v0APts[3+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45;
      v0APts[4+8*i] = (fV0AR4+fV0AFraWd/2.)*sin45+preShapeSec1;
      v0APts[5+8*i] = (fV0AR4+fV0AFraWd/2.)*sin45+2.*preShapeSec1;
      v0APts[6+8*i] = (fV0AR4+fV0AFraWd)*sin45+preShapeSec1;
      v0APts[7+8*i] = fV0AR4*sin45+preShapeSec1;
    }
    new TGeoArb8("sV0ACha2Sec1", fV0ASciWd/2.+2.*preShapeSec1, v0APts);
    new TGeoCompositeShape("sV0ACha12Sec1","sV0ACha1Sec1+sV0ACha2Sec1");
    new TGeoTube("sV0ANail1SciHoleSec1", 0.0, 0.4, 1.65);
    TGeoTranslation *pos1Sec1 = new TGeoTranslation("pos1Sec1", 42.9, 0.51, 0.0);
    pos1Sec1->RegisterYourself();
    new TGeoTube("sV0ANail2SciHoleSec1", 0.0, 0.4, 1.65);
    TGeoTranslation *pos2Sec1 = new TGeoTranslation("pos2Sec1", 30.73,29.98,0.0);
    pos2Sec1->RegisterYourself();
    new TGeoCompositeShape("sV0ANailsSciHolesSec1","sV0ANail1SciHoleSec1:pos1Sec1+sV0ANail2SciHoleSec1:pos2Sec1");
    new TGeoCompositeShape("sV0AChaSec1","sV0ACha12Sec1+sV0ANailsSciHolesSec1");
    new TGeoTubeSeg("sV0AFicR5Sec1", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShapeSec1, 0, 45);
    new TGeoBBox("sV0AFicFEEBSec1", fV0AFEEBWd/2., fV0AFEEBHt/2., fV0AFEEBTh/2.);
    TGeoRotation *rot = new TGeoRotation("rot");
    rot->RotateX(90);
    rot->RotateZ(-90.0);
    double aFEEshiftR2Sec1 = fV0AR6  +  fV0AOctH2 + fV0APlaAl;
    TGeoCombiTrans *posFicFEEBSec1 = new TGeoCombiTrans("posFicFEEBSec1", aFEEshiftR2Sec1*cos225 + 2.0, 0, 7.5, rot);
    posFicFEEBSec1->RegisterYourself();
    rot = new TGeoRotation("rot");
    rot->RotateX(90);
    rot->RotateZ(-90.0+45.0);
    TGeoCombiTrans *posFicFEEBUpSec1 = new TGeoCombiTrans("posFicFEEBUpSec1", (aFEEshiftR2Sec1*cos225 + 2.0 )*cos45, (aFEEshiftR2Sec1*cos225 + 2.0 )*sin45, 7.5, rot);
    posFicFEEBUpSec1->RegisterYourself();
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = fV0AOctH2/2.;	       v0APts[1+8*i] = fV0AFEEBHt/2. + 2.5;
    v0APts[2+8*i] = fV0AOctH2/2.;              v0APts[3+8*i] = -fV0AFEEBHt/2.- 2.5;
    v0APts[4+8*i] = -fV0AOctH2/2.;	       v0APts[5+8*i] = -fV0AFEEBHt/2.- 2.5;
    v0APts[6+8*i] = -fV0AOctH2/2.;	       v0APts[7+8*i] = fV0AFEEBHt/2.+ 2.5;
    }
    new TGeoArb8("sV0AFicOct2Sec1", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoTranslation *posFicOct2Sec1 = new TGeoTranslation("posFicOct2Sec1",(aFEEshiftR2Sec1*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0,0.0,0.0);
    posFicOct2Sec1->RegisterYourself();  
    rot = new TGeoRotation("rot");
    rot->RotateZ(-90.0+45.0+90.0);
    TGeoCombiTrans *posFicOct2UpSec1 = new TGeoCombiTrans("posFicOct2UpSec1",((aFEEshiftR2Sec1*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0)*cos45,((aFEEshiftR2Sec1*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0)*sin45,0.0,rot);
    posFicOct2UpSec1->RegisterYourself(); 

    /// Frame
    TGeoVolume *v0AFraSec1 = new TGeoVolumeAssembly("V0AFraSec1");
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[3+8*i] = fV0AFraWd/2.;
      v0APts[4+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[5+8*i] = fV0AFraWd/2.;
      v0APts[6+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[7+8*i] = 0.;
    }
    TGeoArb8 *sV0AFraB1Sec1 = new TGeoArb8("sV0AFraB1Sec1",fV0ASciWd/2.,v0APts);
    TGeoVolume *v0AFraB1Sec1 = new TGeoVolume("V0AFraB1Sec1",sV0AFraB1Sec1,medV0AFra);
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
    TGeoArb8 *sV0AFraB2Sec1 = new TGeoArb8("sV0AFraB2Sec1", fV0ASciWd/2., v0APts);
    TGeoVolume *v0AFraB2Sec1 = new TGeoVolume("V0AFraB2Sec1",sV0AFraB2Sec1,medV0AFra);
    v0AFraB1Sec1->SetLineColor(kV0AColorFra); v0AFraB2Sec1->SetLineColor(kV0AColorFra);
    v0AFraSec1->AddNode(v0AFraB1Sec1,1);
    v0AFraSec1->AddNode(v0AFraB2Sec1,1);  
    new TGeoTubeSeg( "sV0AFraR1bSec1", fV0AR0-fV0AFraWd/2.,
		     fV0AR0+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR2bSec1", fV0AR1-fV0AFraWd/2.,
		     fV0AR1+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR3bSec1", fV0AR2-fV0AFraWd/2.,
		     fV0AR2+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR4bSec1", fV0AR3-fV0AFraWd/2.,
		     fV0AR3+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR5bSec1", fV0AR4-fV0AFraWd/2.,
		     fV0AR4+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    TGeoCompositeShape *sV0AFraR1Sec1 = new TGeoCompositeShape("sV0AFraR1Sec1","sV0AFraR1bSec1-sV0AChaSec1");
    TGeoCompositeShape *sV0AFraR2Sec1 = new TGeoCompositeShape("sV0AFraR2Sec1","sV0AFraR2bSec1-sV0AChaSec1");
    TGeoCompositeShape *sV0AFraR3Sec1 = new TGeoCompositeShape("sV0AFraR3Sec1","sV0AFraR3bSec1-sV0AChaSec1");
    TGeoCompositeShape *sV0AFraR4Sec1 = new TGeoCompositeShape("sV0AFraR4Sec1","sV0AFraR4bSec1-sV0AChaSec1");
    TGeoCompositeShape *sV0AFraR5Sec1 = new TGeoCompositeShape("sV0AFraR5Sec1","sV0AFraR5bSec1-sV0AChaSec1");
    TGeoVolume *v0AFraR1Sec1 = new TGeoVolume("V0AFraR1Sec1",sV0AFraR1Sec1,medV0AFra);
    TGeoVolume *v0AFraR2Sec1 = new TGeoVolume("V0AFraR2Sec1",sV0AFraR2Sec1,medV0AFra);
    TGeoVolume *v0AFraR3Sec1 = new TGeoVolume("V0AFraR3Sec1",sV0AFraR3Sec1,medV0AFra);
    TGeoVolume *v0AFraR4Sec1 = new TGeoVolume("V0AFraR4Sec1",sV0AFraR4Sec1,medV0AFra);
    TGeoVolume *v0AFraR5Sec1 = new TGeoVolume("V0AFraR5Sec1",sV0AFraR5Sec1,medV0AFra);
    v0AFraR1Sec1->SetLineColor(kV0AColorFra); v0AFraR2Sec1->SetLineColor(kV0AColorFra);
    v0AFraR3Sec1->SetLineColor(kV0AColorFra); v0AFraR4Sec1->SetLineColor(kV0AColorFra);
    v0AFraR5Sec1->SetLineColor(kV0AColorFra);
    v0AFraSec1->AddNode(v0AFraR1Sec1,1); 
    v0AFraSec1->AddNode(v0AFraR2Sec1,1);
    v0AFraSec1->AddNode(v0AFraR3Sec1,1); 
    v0AFraSec1->AddNode(v0AFraR4Sec1,1);
    v0AFraSec1->AddNode(v0AFraR5Sec1,1);
    v0ASec1->AddNode(v0AFraSec1,1);
    
    /// Sensitive scintilator
    TGeoVolume *v0ASciSec1 = new TGeoVolumeAssembly("V0ASciSec1");
    new TGeoTubeSeg( "sV0AR1bSec1", fV0AR0+fV0AFraWd/2.,
		     fV0AR1-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR2bSec1", fV0AR1+fV0AFraWd/2.,
		     fV0AR2-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR3bSec1", fV0AR2+fV0AFraWd/2.,
		     fV0AR3-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR4bSec1", fV0AR3+fV0AFraWd/2.,
		     fV0AR4-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    TGeoCompositeShape *sV0AR1Sec1 = new TGeoCompositeShape("sV0AR1Sec1","sV0AR1bSec1-sV0AChaSec1");
    TGeoCompositeShape *sV0AR2Sec1 = new TGeoCompositeShape("sV0AR2Sec1","sV0AR2bSec1-sV0AChaSec1");
    TGeoCompositeShape *sV0AR3Sec1 = new TGeoCompositeShape("sV0AR3Sec1","sV0AR3bSec1-sV0AChaSec1");
    TGeoCompositeShape *sV0AR4Sec1 = new TGeoCompositeShape("sV0AR4Sec1","sV0AR4bSec1-sV0AChaSec1");
    TGeoVolume *v0L1Sec1 = new TGeoVolume("V0L1Sec1",sV0AR1Sec1,medV0ASci);
    TGeoVolume *v0L2Sec1 = new TGeoVolume("V0L2Sec1",sV0AR2Sec1,medV0ASci);
    TGeoVolume *v0L3Sec1 = new TGeoVolume("V0L3Sec1",sV0AR3Sec1,medV0ASci);
    TGeoVolume *v0L4Sec1 = new TGeoVolume("V0L4Sec1",sV0AR4Sec1,medV0ASci);
    v0L1Sec1->SetLineColor(kV0AColorSci); v0L2Sec1->SetLineColor(kV0AColorSci);
    v0L3Sec1->SetLineColor(kV0AColorSci); v0L4Sec1->SetLineColor(kV0AColorSci);
    v0ASec1->AddNode(v0L1Sec1,1);
    v0ASec1->AddNode(v0L2Sec1,1);
    v0ASec1->AddNode(v0L3Sec1,1);
    v0ASec1->AddNode(v0L4Sec1,1);      
    
    /// Segment of octagon 
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] =  fV0AR6-fV0AOctH2;           v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = (fV0AR6-fV0AOctH2)*sin45;  v0APts[3+8*i] = (fV0AR6-fV0AOctH2)*sin45;
      v0APts[4+8*i] = fV0AR6*sin45;		 v0APts[5+8*i] = fV0AR6*sin45;
      v0APts[6+8*i] = fV0AR6;			 v0APts[7+8*i] = 0.;
    }
    new TGeoArb8("sV0AOct2Sec1", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoCompositeShape *sV0AOct2FEEBSec1 = new TGeoCompositeShape("sV0AOct2FEEBSec1","sV0AOct2Sec1-sV0AFicFEEBSec1:posFicFEEBSec1-sV0AFicFEEBSec1:posFicFEEBUpSec1-sV0AFicOct2Sec1:posFicOct2Sec1-sV0AFicOct2Sec1:posFicOct2UpSec1");
    TGeoVolume *v0AOct2Sec1 = new TGeoVolume("V0AOct2Sec1", sV0AOct2FEEBSec1,medV0ASup);
    v0AOct2Sec1->SetLineColor(kV0AColorOct);
    TGeoVolume *v0ASupSec1 = new TGeoVolumeAssembly("V0ASupSec1");
    v0ASupSec1->AddNode(v0AOct2Sec1,1);
    v0ASec1->AddNode(v0ASupSec1,1);

    //Bunch of fibers
    v0APts[ 0] = v0APts[ 2] = -13.0;
    v0APts[ 1] = v0APts[ 7] = (fV0ASciWd+fV0AOctWd)/2.-0.01;
    v0APts[ 3] = v0APts[ 5] = (fV0ASciWd+fV0AOctWd)/2.+0.01;
    v0APts[ 4] = v0APts[ 6] = +13.0;
    v0APts[ 8] = v0APts[10] = -10.0;
    v0APts[ 9] = v0APts[15] = 0.;
    v0APts[11] = v0APts[13] = 0.25;
    v0APts[12] = v0APts[14] = +10.0;
    new TGeoArb8("sV0AFib1Sec1", 11.5, v0APts);
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateZ(-90.+22.5);
    TGeoCombiTrans *fibpos1Sec1 = new TGeoCombiTrans("fibpos1Sec1", (fV0AR6-fV0AOctH2+fV0AR5)*cos225/2. - 3.8, (fV0AR6-fV0AOctH2+fV0AR5)*sin225/2.-1.8, 0, rot);
    fibpos1Sec1->RegisterYourself();
    TGeoCompositeShape *sV0AFib1HoleSec1 = new TGeoCompositeShape("sV0AFib1HoleSec1","sV0AFib1Sec1:fibpos1Sec1-sV0AFicR5Sec1"); 
    TGeoVolume *v0AFib1HoleSec1 = new TGeoVolume("V0AFib1HoleSec1",sV0AFib1HoleSec1,medV0AFib);
    v0AFib1HoleSec1->SetLineColor(kV0AColorFib);
    new TGeoArb8("sV0AFib2Sec1", 11.5, v0APts);
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateY(180);
    rot->RotateZ(-90.+22.5);
    TGeoCombiTrans *fibpos2Sec1 = new TGeoCombiTrans("fibpos2Sec1", (fV0AR6-fV0AOctH2+fV0AR5)*cos225/2. - 3.8, (fV0AR6-fV0AOctH2+fV0AR5)*sin225/2.-1.8, 0, rot);
    fibpos2Sec1->RegisterYourself();
    TGeoCompositeShape *sV0AFib2HoleSec1 = new TGeoCompositeShape("sV0AFib2HoleSec1","sV0AFib2Sec1:fibpos2Sec1-sV0AFicR5Sec1");
    TGeoVolume *v0AFib2HoleSec1 = new TGeoVolume("V0AFib2HoleSec1",sV0AFib2HoleSec1,medV0AFib);
    v0AFib2HoleSec1->SetLineColor(kV0AColorFib);
    TGeoVolume *v0AFibSec1 = new TGeoVolumeAssembly("V0AFibSec1");
    v0AFibSec1->AddNode(v0AFib1HoleSec1,1);
    v0AFibSec1->AddNode(v0AFib2HoleSec1,1);
    v0ASec1->AddNode(v0AFibSec1,1); 
    
     /// Plates
    new TGeoTube("sV0ANail1PlaInHoleSec1", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoTube("sV0ANail2PlaInHoleSec1", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoCompositeShape("sV0ANailsPlaInHolesSec1","sV0ANail1PlaInHoleSec1:pos1Sec1+sV0ANail2PlaInHoleSec1:pos2Sec1");
    new TGeoTube("sV0ANail1PlaOuHoleSec1", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoTube("sV0ANail2PlaOuHoleSec1", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoCompositeShape("sV0ANailsPlaOuHolesSec1","sV0ANail1PlaOuHoleSec1:pos1Sec1+sV0ANail2PlaOuHoleSec1:pos2Sec1");
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0;			v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = fV0AR0*sin45;		v0APts[3+8*i] = fV0AR0*sin45;
      v0APts[4+8*i] = fV0AR6 * sin45;	v0APts[5+8*i] = fV0AR6*sin45;
      v0APts[6+8*i] = fV0AR6;		v0APts[7+8*i] = 0.;
    }
    new TGeoArb8("sV0APlaInSec1", (fV0APlaWd-2*fV0APlaAl)/2., v0APts);
    TGeoCompositeShape *sV0APlaInNailsHolesSec1 = new TGeoCompositeShape("sV0APlaInNailsHolesSec1","sV0APlaInSec1-sV0ANailsPlaInHolesSec1-sV0AFicFEEBSec1:posFicFEEBSec1-sV0AFicFEEBSec1:posFicFEEBUpSec1");
    TGeoVolume *v0APlaInNailsHolesSec1 = new TGeoVolume("V0APlaInNailsHolesSec1", sV0APlaInNailsHolesSec1, medV0APlaIn);
    new TGeoArb8("sV0APlaOuSec1", fV0APlaAl/2., v0APts);
    TGeoCompositeShape *sV0APlaOuNailsHolesSec1 = new TGeoCompositeShape("sV0APlaOuNailsHolesSec1","sV0APlaOuSec1-sV0ANailsPlaOuHolesSec1-sV0AFicFEEBSec1:posFicFEEBSec1-sV0AFicFEEBSec1:posFicFEEBUpSec1"); 
    TGeoVolume *v0APlaOuNailsHolesSec1 = new TGeoVolume("V0APlaOuNailsHolesSec1", sV0APlaOuNailsHolesSec1, medV0APlaOu);
    v0APlaInNailsHolesSec1->SetLineColor(kV0AColorPlaIn); v0APlaOuNailsHolesSec1->SetLineColor(kV0AColorPlaOu);
    TGeoVolume *v0APlaSec1 = new TGeoVolumeAssembly("V0APlaSec1");
    v0APlaSec1->AddNode(v0APlaInNailsHolesSec1,1);
    v0APlaSec1->AddNode(v0APlaOuNailsHolesSec1,1,new TGeoTranslation(0,0,(fV0APlaWd-fV0APlaAl)/2.));
    v0APlaSec1->AddNode(v0APlaOuNailsHolesSec1,2,new TGeoTranslation(0,0,-(fV0APlaWd-fV0APlaAl)/2.));
    v0ASec1->AddNode(v0APlaSec1,1,new TGeoTranslation(0,0,(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    v0ASec1->AddNode(v0APlaSec1,2,new TGeoTranslation(0,0,-(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    
     /// Non-sensitive scintilator
    new TGeoTubeSeg("sV0AR5S2Sec1", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShapeSec1, 0, 45);
    TGeoCompositeShape *sV0AR5Sec1 = new TGeoCompositeShape("V0AR5Sec1","sV0AR5S2Sec1 - sV0AChaSec1");
    TGeoVolume *v0AR5Sec1 = new TGeoVolume("V0AR5Sec1",sV0AR5Sec1,medV0ASci);
    v0AR5Sec1->SetLineColor(kV0AColorSci);
    v0ASciSec1->AddNode(v0AR5Sec1,1);
    v0ASec1->AddNode(v0ASciSec1,1); 

    /// PMBox
    TGeoVolume* v0APMSec1 = new TGeoVolumeAssembly("V0APMSec1");
    new TGeoBBox("sV0APMB1Sec1", fV0APMBWd/2., fV0APMBHt/2., fV0APMBTh/2.);
    new TGeoBBox("sV0APMB2Sec1", fV0APMBWd/2.-fV0APMBWdW, fV0APMBHt/2.-fV0APMBHtW, fV0APMBTh/2.-fV0APMBThW);
    TGeoCompositeShape *sV0APMBSec1 = new TGeoCompositeShape("sV0APMBSec1","sV0APMB1Sec1-sV0APMB2Sec1");
    TGeoVolume *v0APMBSec1 = new TGeoVolume("V0APMBSec1",sV0APMBSec1, medV0APMAlum);
    v0APMBSec1->SetLineColor(kV0AColorPMA);
    v0APMSec1->AddNode(v0APMBSec1,1);

    /// PMTubes
    TGeoTube *sV0APMT1Sec1 = new TGeoTube("sV0APMT1Sec1", fV0APMTR1, fV0APMTR2, fV0APMTH/2.);
    TGeoVolume *v0APMT1Sec1 = new TGeoVolume("V0APMT1Sec1", sV0APMT1Sec1, medV0APMGlass);
    TGeoTube *sV0APMT2Sec1 = new TGeoTube("sV0APMT2Sec1", fV0APMTR3, fV0APMTR4, fV0APMTH/2.);
    TGeoVolume *v0APMT2Sec1 = new TGeoVolume("V0APMT2Sec1", sV0APMT2Sec1, medV0APMAlum);
    TGeoVolume *v0APMTSec1 = new TGeoVolumeAssembly("V0APMTSec1");
    TGeoTube *sV0APMTTSec1 = new TGeoTube("sV0APMTTSec1", 0., fV0APMTR4, fV0APMTB/2.);
    TGeoVolume *v0APMTTSec1 = new TGeoVolume("V0APMTTSec1", sV0APMTTSec1, medV0APMAlum);
    v0APMT1Sec1->SetLineColor(kV0AColorPMG);
    v0APMT2Sec1->SetLineColor(kV0AColorPMA);
    v0APMTTSec1->SetLineColor(kV0AColorPMA);
    rot = new TGeoRotation("rot", 90, 0, 180, 0, 90, 90);
    v0APMTSec1->AddNode(v0APMT1Sec1,1,rot);
    v0APMTSec1->AddNode(v0APMT2Sec1,1,rot);
    v0APMTSec1->AddNode(v0APMTTSec1,1,new TGeoCombiTrans(0,-(fV0APMTH+fV0APMTB)/2.,0,rot));
    double autoShiftSec1 = (fV0APMBWd-2*fV0APMBWdW)/4.;
    v0APMSec1->AddNode(v0APMTSec1, 1, new TGeoTranslation(-1.5*autoShiftSec1, 0, 0));
    v0APMSec1->AddNode(v0APMTSec1, 2, new TGeoTranslation(-0.5*autoShiftSec1, 0, 0));
    v0APMSec1->AddNode(v0APMTSec1, 3, new TGeoTranslation(+0.5*autoShiftSec1, 0, 0));
    v0APMSec1->AddNode(v0APMTSec1, 4, new TGeoTranslation(+1.5*autoShiftSec1, 0, 0));

    // PM
    rot = new TGeoRotation("rot");
    rot->RotateX(90-fV0APMBAng);
    rot->RotateZ(-90.+22.5);
    double cosAngPMBSec1 = TMath::Cos(fV0APMBAng*TMath::DegToRad());
    double sinAngPMBSec1 = TMath::Sin(fV0APMBAng*TMath::DegToRad());
    double shiftZSec1 = fV0APMBHt/2. * cosAngPMBSec1
      -   ( fV0ASciWd + 2 * fV0AOctWd + 2 * fV0APlaWd )/2.   -   fV0APMBTh/2. * sinAngPMBSec1;
    double shiftRSec1 = fV0AR6  +  fV0AOctH2 + fV0APlaAl;
    v0ASec1->AddNode(v0APMSec1,1, new TGeoCombiTrans( shiftRSec1*cos225+1.07, shiftRSec1*sin225, shiftZSec1, rot));
    
    // Aluminium nails 
    TGeoTube *sV0ANail1Sec1 = new TGeoTube("sV0ANail1Sec1", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail1Sec1 = new TGeoVolume("V0ANail1Sec1", sV0ANail1Sec1, medV0APMAlum);
    v0ANail1Sec1->SetLineColor(kV0AColorPMA);// this is the color for aluminium
    v0ASec1->AddNode(v0ANail1Sec1,1,new TGeoTranslation(42.9, 0.51, 0.0));
    TGeoTube *sV0ANail2Sec1 = new TGeoTube("sV0ANail2Sec1", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail2Sec1 = new TGeoVolume("V0ANail2Sec1", sV0ANail2Sec1, medV0APMAlum);
    v0ANail2Sec1->SetLineColor(kV0AColorPMA);
    v0ASec1->AddNode(v0ANail2Sec1,1,new TGeoTranslation(30.73,29.98,0.0)); 
        
    /// Adding sector to v0LE volume
    for(int i=0; i<1; i++) {
       TGeoRotation *rotation = new TGeoRotation("rotation", 90., i*45., 90., 90.+i*45., 0., 0.);
       v0LE->AddNode(v0ASec1,i+1,rotation);  
    }

    //Front end electronics for sector 1

    //FEEBox
    TGeoVolume* v0AFEE = new TGeoVolumeAssembly("V0AFEE");
    new TGeoBBox("sV0AFEEB1", fV0AFEEBWd/2., fV0AFEEBHt/2., fV0AFEEBTh/2.);
    new TGeoBBox("sV0AFEEB2", fV0AFEEBWd/2.-fV0APMBWdW, fV0AFEEBHt/2.-fV0APMBHtW, fV0AFEEBTh/2.-fV0APMBThW);
    TGeoCompositeShape *sV0AFEEB = new TGeoCompositeShape("sV0AFEEB","sV0AFEEB1-sV0AFEEB2");
    TGeoVolume *v0AFEEB = new TGeoVolume("V0AFEEB",sV0AFEEB, medV0APMAlum);
    v0AFEEB->SetLineColor(kV0AColorPMA);
    v0AFEE->AddNode(v0AFEEB,1);

    //Mother and daughter boards
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = -fV0APMBHtW/10.;	       v0APts[1+8*i] = fV0AFEEBTh/2.-fV0APMTB;
    v0APts[2+8*i] = fV0APMBHtW/10.;             v0APts[3+8*i] = fV0AFEEBTh/2.-fV0APMTB;
    v0APts[4+8*i] = fV0APMBHtW/10.;	       v0APts[5+8*i] = -fV0AFEEBTh/2.+fV0APMTB;
    v0APts[6+8*i] = -fV0APMBHtW/10.;	       v0APts[7+8*i] = -fV0AFEEBTh/2.+fV0APMTB;
    }
    TGeoArb8 *sV0AFEEDaughter = new TGeoArb8("sV0AFEEDaughter", fV0AFEEBTh/2.-fV0APMTB, v0APts);
    TGeoVolume *v0AFEEDaughter = new TGeoVolume("V0AFEEDaughter", sV0AFEEDaughter, medV0AFibGlass);
    v0AFEEDaughter->SetLineColor(kV0AColorFibGlass);
    double spacing = fV0APMBHtW;
    v0AFEE->AddNode(v0AFEEDaughter, 1, new TGeoTranslation(9.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE->AddNode(v0AFEEDaughter, 2, new TGeoTranslation(6.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE->AddNode(v0AFEEDaughter, 3, new TGeoTranslation(3.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE->AddNode(v0AFEEDaughter, 4, new TGeoTranslation(0.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE->AddNode(v0AFEEDaughter, 5, new TGeoTranslation(-3.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE->AddNode(v0AFEEDaughter, 6, new TGeoTranslation(-6.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE->AddNode(v0AFEEDaughter, 7, new TGeoTranslation(-9.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE->AddNode(v0AFEEDaughter, 8, new TGeoTranslation(-12.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));   
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = -fV0AFEEBWd/2.+fV0APMBWdW;	       v0APts[1+8*i] = fV0AFEEBTh+fV0APMTB;
    v0APts[2+8*i] = fV0AFEEBWd/2.-fV0APMBWdW;          v0APts[3+8*i] = fV0AFEEBTh+fV0APMTB;
    v0APts[4+8*i] = fV0AFEEBWd/2.-fV0APMBWdW;	       v0APts[5+8*i] = -fV0AFEEBTh-fV0APMTB;
    v0APts[6+8*i] = -fV0AFEEBWd/2.+fV0APMBWdW;	       v0APts[7+8*i] = -fV0AFEEBTh-fV0APMTB;
    }
    TGeoArb8 *sV0AFEEMother = new TGeoArb8("sV0AFEEMother", fV0APMBHtW/10., v0APts);
    TGeoVolume *v0AFEEMother = new TGeoVolume("V0AFEEMother", sV0AFEEMother, medV0AFibGlass);
    v0AFEEMother->SetLineColor(kV0AColorFibGlass);
    v0AFEE->AddNode(v0AFEEMother, 1, new TGeoTranslation(0.0, 0.0, -fV0AFEEBTh/2.+fV0APMBThW+fV0APMBHtW));
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = -fV0AFEEBWd/2.+fV0APMBWdW;	       v0APts[1+8*i] = (fV0AFEEBTh+fV0APMTB)/2.;
      v0APts[2+8*i] = fV0AFEEBWd/2.-fV0APMBWdW;                v0APts[3+8*i] = (fV0AFEEBTh+fV0APMTB)/2.;
      v0APts[4+8*i] = fV0AFEEBWd/2.-fV0APMBWdW;	               v0APts[5+8*i] = (-fV0AFEEBTh-fV0APMTB)/2.;
      v0APts[6+8*i] = -fV0AFEEBWd/2.+fV0APMBWdW;	       v0APts[7+8*i] = (-fV0AFEEBTh-fV0APMTB)/2.; 
    }
    TGeoArb8 *sV0AFEEHalfMother = new TGeoArb8("sV0AFEEHalfMother", fV0APMBHtW/10., v0APts);
    TGeoVolume *v0AFEEHalfMother = new TGeoVolume("V0AFEEHalfMother", sV0AFEEHalfMother, medV0AFibGlass);
    v0AFEEHalfMother->SetLineColor(kV0AColorFibGlass);
    v0AFEE->AddNode(v0AFEEHalfMother, 1, new TGeoTranslation(0.0, -(fV0AFEEBTh+fV0APMTB)/2., -2.0*spacing));


    //FEE
    rot = new TGeoRotation("rot");
    rot->RotateX(90);
    rot->RotateZ(+90.0);
    double aFEEshiftR = fV0AR6  +  fV0AOctH2 + fV0APlaAl;
    v0LE->AddNode(v0AFEE,1, new TGeoCombiTrans( aFEEshiftR*cos225+2.0, 0, 7.5, rot));
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = fV0AOctH2/2.;	       v0APts[1+8*i] = fV0AFEEBHt/2. + 2.5;
    v0APts[2+8*i] = fV0AOctH2/2.;              v0APts[3+8*i] = -fV0AFEEBHt/2.- 2.5;
    v0APts[4+8*i] = -fV0AOctH2/2.;	       v0APts[5+8*i] = -fV0AFEEBHt/2.- 2.5;
    v0APts[6+8*i] = -fV0AOctH2/2.;	       v0APts[7+8*i] = fV0AFEEBHt/2.+ 2.5;
    }
    TGeoArb8 *sV0AFEEOct2 = new TGeoArb8("sV0AFEEOct2", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoVolume *v0AFEEOct2 = new TGeoVolume("V0AFEEOct2",sV0AFEEOct2, medV0ASup);
    v0AFEEOct2->SetLineColor(kV0AColorOct);
    v0LE->AddNode(v0AFEEOct2,1, new TGeoTranslation((aFEEshiftR2Sec1*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0,0.0,0.0));


    /// Definition sector 2
    TGeoVolume *v0ASec2 = new TGeoVolumeAssembly("V0ASec2");
        
    /// For boolean sustraction
    double preShapeSec2 = 0.2;
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0-fV0AFraWd/2.-preShapeSec2;  v0APts[1+8*i] = -preShapeSec2;
      v0APts[2+8*i] = fV0AR0-fV0AFraWd/2.-preShapeSec2;  v0APts[3+8*i] = fV0AFraWd/2.;
      v0APts[4+8*i] = fV0AR4+fV0AFraWd/2.+preShapeSec2;  v0APts[5+8*i] = fV0AFraWd/2.;
      v0APts[6+8*i] = fV0AR4+fV0AFraWd/2.+preShapeSec2;  v0APts[7+8*i] = -preShapeSec2;
    }
    new TGeoArb8("sV0ACha1Sec2",fV0ASciWd/1.5,v0APts);
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0*sin45-preShapeSec2;
      v0APts[1+8*i] = (fV0AR0-fV0AFraWd)*sin45-preShapeSec2;
      v0APts[2+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45-preShapeSec2;
      v0APts[3+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45;
      v0APts[4+8*i] = (fV0AR4+fV0AFraWd/2.)*sin45+preShapeSec2;
      v0APts[5+8*i] = (fV0AR4+fV0AFraWd/2.)*sin45+2.*preShapeSec2;
      v0APts[6+8*i] = (fV0AR4+fV0AFraWd)*sin45+preShapeSec2;
      v0APts[7+8*i] = fV0AR4*sin45+preShapeSec2;
    }
    new TGeoArb8("sV0ACha2Sec2", fV0ASciWd/2.+2.*preShapeSec2, v0APts);
    new TGeoCompositeShape("sV0ACha12Sec2","sV0ACha1Sec2+sV0ACha2Sec2");
    new TGeoTube("sV0ANail1SciHoleSec2", 0.0, 0.4, 1.65);
    TGeoTranslation *pos1Sec2 = new TGeoTranslation("pos1Sec2", 42.9, 0.51, 0.0);
    pos1Sec2->RegisterYourself();
    new TGeoTube("sV0ANail2SciHoleSec2", 0.0, 0.4, 1.65);
    TGeoTranslation *pos2Sec2 = new TGeoTranslation("pos2Sec2", 30.73,29.98,0.0);
    pos2Sec2->RegisterYourself();
    new TGeoCompositeShape("sV0ANailsSciHolesSec2","sV0ANail1SciHoleSec2:pos1Sec2+sV0ANail2SciHoleSec2:pos2Sec2");
    new TGeoCompositeShape("sV0AChaSec2","sV0ACha12Sec2+sV0ANailsSciHolesSec2");
    new TGeoTubeSeg("sV0AFicR5Sec2", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShapeSec2, 0, 45);
    new TGeoBBox("sV0AFicFEEBSec2", fV0AFEEBWd/2., fV0AFEEBHt/2., fV0AFEEBTh/2.);
    rot = new TGeoRotation("rot");
    rot->RotateX(90);
    rot->RotateZ(-90.0);
    double aFEEshiftR2Sec2 = fV0AR6  +  fV0AOctH2 + fV0APlaAl;
    TGeoCombiTrans *posFicFEEBSec2 = new TGeoCombiTrans("posFicFEEBSec2", aFEEshiftR2Sec2*cos225 + 2.0, 0, 7.5, rot);
    posFicFEEBSec2->RegisterYourself();
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = fV0AOctH2/2.;	       v0APts[1+8*i] = fV0AFEEBHt/2. + 2.5;
    v0APts[2+8*i] = fV0AOctH2/2.;              v0APts[3+8*i] = -fV0AFEEBHt/2.- 2.5;
    v0APts[4+8*i] = -fV0AOctH2/2.;	       v0APts[5+8*i] = -fV0AFEEBHt/2.- 2.5;
    v0APts[6+8*i] = -fV0AOctH2/2.;	       v0APts[7+8*i] = fV0AFEEBHt/2.+ 2.5;
    }
    new TGeoArb8("sV0AFicOct2Sec2", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoTranslation *posFicOct2Sec2 = new TGeoTranslation("posFicOct2Sec2",(aFEEshiftR2Sec2*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0,0.0,0.0);
    posFicOct2Sec2->RegisterYourself();  

    /// Frame
    TGeoVolume *v0AFraSec2 = new TGeoVolumeAssembly("V0AFraSec2");
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[3+8*i] = fV0AFraWd/2.;
      v0APts[4+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[5+8*i] = fV0AFraWd/2.;
      v0APts[6+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[7+8*i] = 0.;
    }
    TGeoArb8 *sV0AFraB1Sec2 = new TGeoArb8("sV0AFraB1Sec2",fV0ASciWd/2.,v0APts);
    TGeoVolume *v0AFraB1Sec2 = new TGeoVolume("V0AFraB1Sec2",sV0AFraB1Sec2,medV0AFra);
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
    TGeoArb8 *sV0AFraB2Sec2 = new TGeoArb8("sV0AFraB2Sec2", fV0ASciWd/2., v0APts);
    TGeoVolume *v0AFraB2Sec2 = new TGeoVolume("V0AFraB2Sec2",sV0AFraB2Sec2,medV0AFra);
    v0AFraB1Sec2->SetLineColor(kV0AColorFra); v0AFraB2Sec2->SetLineColor(kV0AColorFra);
    v0AFraSec2->AddNode(v0AFraB1Sec2,1);
    v0AFraSec2->AddNode(v0AFraB2Sec2,1);  
    new TGeoTubeSeg( "sV0AFraR1bSec2", fV0AR0-fV0AFraWd/2.,
		     fV0AR0+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR2bSec2", fV0AR1-fV0AFraWd/2.,
		     fV0AR1+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR3bSec2", fV0AR2-fV0AFraWd/2.,
		     fV0AR2+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR4bSec2", fV0AR3-fV0AFraWd/2.,
		     fV0AR3+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR5bSec2", fV0AR4-fV0AFraWd/2.,
		     fV0AR4+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    TGeoCompositeShape *sV0AFraR1Sec2 = new TGeoCompositeShape("sV0AFraR1Sec2","sV0AFraR1bSec2-sV0AChaSec2");
    TGeoCompositeShape *sV0AFraR2Sec2 = new TGeoCompositeShape("sV0AFraR2Sec2","sV0AFraR2bSec2-sV0AChaSec2");
    TGeoCompositeShape *sV0AFraR3Sec2 = new TGeoCompositeShape("sV0AFraR3Sec2","sV0AFraR3bSec2-sV0AChaSec2");
    TGeoCompositeShape *sV0AFraR4Sec2 = new TGeoCompositeShape("sV0AFraR4Sec2","sV0AFraR4bSec2-sV0AChaSec2");
    TGeoCompositeShape *sV0AFraR5Sec2 = new TGeoCompositeShape("sV0AFraR5Sec2","sV0AFraR5bSec2-sV0AChaSec2");
    TGeoVolume *v0AFraR1Sec2 = new TGeoVolume("V0AFraR1Sec2",sV0AFraR1Sec2,medV0AFra);
    TGeoVolume *v0AFraR2Sec2 = new TGeoVolume("V0AFraR2Sec2",sV0AFraR2Sec2,medV0AFra);
    TGeoVolume *v0AFraR3Sec2 = new TGeoVolume("V0AFraR3Sec2",sV0AFraR3Sec2,medV0AFra);
    TGeoVolume *v0AFraR4Sec2 = new TGeoVolume("V0AFraR4Sec2",sV0AFraR4Sec2,medV0AFra);
    TGeoVolume *v0AFraR5Sec2 = new TGeoVolume("V0AFraR5Sec2",sV0AFraR5Sec2,medV0AFra);
    v0AFraR1Sec2->SetLineColor(kV0AColorFra); v0AFraR2Sec2->SetLineColor(kV0AColorFra);
    v0AFraR3Sec2->SetLineColor(kV0AColorFra); v0AFraR4Sec2->SetLineColor(kV0AColorFra);
    v0AFraR5Sec2->SetLineColor(kV0AColorFra);
    v0AFraSec2->AddNode(v0AFraR1Sec2,1); 
    v0AFraSec2->AddNode(v0AFraR2Sec2,1);
    v0AFraSec2->AddNode(v0AFraR3Sec2,1); 
    v0AFraSec2->AddNode(v0AFraR4Sec2,1);
    v0AFraSec2->AddNode(v0AFraR5Sec2,1);
    v0ASec2->AddNode(v0AFraSec2,1);
    
    /// Sensitive scintilator
    TGeoVolume *v0ASciSec2 = new TGeoVolumeAssembly("V0ASciSec2");
    new TGeoTubeSeg( "sV0AR1bSec2", fV0AR0+fV0AFraWd/2.,
		     fV0AR1-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR2bSec2", fV0AR1+fV0AFraWd/2.,
		     fV0AR2-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR3bSec2", fV0AR2+fV0AFraWd/2.,
		     fV0AR3-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR4bSec2", fV0AR3+fV0AFraWd/2.,
		     fV0AR4-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    TGeoCompositeShape *sV0AR1Sec2 = new TGeoCompositeShape("sV0AR1Sec2","sV0AR1bSec2-sV0AChaSec2");
    TGeoCompositeShape *sV0AR2Sec2 = new TGeoCompositeShape("sV0AR2Sec2","sV0AR2bSec2-sV0AChaSec2");
    TGeoCompositeShape *sV0AR3Sec2 = new TGeoCompositeShape("sV0AR3Sec2","sV0AR3bSec2-sV0AChaSec2");
    TGeoCompositeShape *sV0AR4Sec2 = new TGeoCompositeShape("sV0AR4Sec2","sV0AR4bSec2-sV0AChaSec2");
    TGeoVolume *v0L1Sec2 = new TGeoVolume("V0L1Sec2",sV0AR1Sec2,medV0ASci);
    TGeoVolume *v0L2Sec2 = new TGeoVolume("V0L2Sec2",sV0AR2Sec2,medV0ASci);
    TGeoVolume *v0L3Sec2 = new TGeoVolume("V0L3Sec2",sV0AR3Sec2,medV0ASci);
    TGeoVolume *v0L4Sec2 = new TGeoVolume("V0L4Sec2",sV0AR4Sec2,medV0ASci);
    v0L1Sec2->SetLineColor(kV0AColorSci); v0L2Sec2->SetLineColor(kV0AColorSci);
    v0L3Sec2->SetLineColor(kV0AColorSci); v0L4Sec2->SetLineColor(kV0AColorSci);
    v0ASec2->AddNode(v0L1Sec2,1);
    v0ASec2->AddNode(v0L2Sec2,1);
    v0ASec2->AddNode(v0L3Sec2,1);
    v0ASec2->AddNode(v0L4Sec2,1);      
    
    /// Segment of octagon 
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] =  fV0AR6-fV0AOctH2;           v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = (fV0AR6-fV0AOctH2)*sin45;  v0APts[3+8*i] = (fV0AR6-fV0AOctH2)*sin45;
      v0APts[4+8*i] = fV0AR6*sin45;		 v0APts[5+8*i] = fV0AR6*sin45;
      v0APts[6+8*i] = fV0AR6;			 v0APts[7+8*i] = 0.;
    }
    new TGeoArb8("sV0AOct2Sec2", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoCompositeShape *sV0AOct2FEEBSec2 = new TGeoCompositeShape("sV0AOct2FEEBSec2","sV0AOct2Sec2-sV0AFicFEEBSec2:posFicFEEBSec2-sV0AFicOct2Sec2:posFicOct2Sec2");
    TGeoVolume *v0AOct2Sec2 = new TGeoVolume("V0AOct2Sec2", sV0AOct2FEEBSec2,medV0ASup);
    v0AOct2Sec2->SetLineColor(kV0AColorOct);
    TGeoVolume *v0ASupSec2 = new TGeoVolumeAssembly("V0ASupSec2");
    v0ASupSec2->AddNode(v0AOct2Sec2,1);
    v0ASec2->AddNode(v0ASupSec2,1);

    //Bunch of fibers
    v0APts[ 0] = v0APts[ 2] = -13.0;
    v0APts[ 1] = v0APts[ 7] = (fV0ASciWd+fV0AOctWd)/2.-0.01;
    v0APts[ 3] = v0APts[ 5] = (fV0ASciWd+fV0AOctWd)/2.+0.01;
    v0APts[ 4] = v0APts[ 6] = +13.0;
    v0APts[ 8] = v0APts[10] = -10.0;
    v0APts[ 9] = v0APts[15] = 0.;
    v0APts[11] = v0APts[13] = 0.25;
    v0APts[12] = v0APts[14] = +10.0;
    new TGeoArb8("sV0AFib1Sec2", 11.5, v0APts);
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateZ(-90.+22.5);
    TGeoCombiTrans *fibpos1Sec2 = new TGeoCombiTrans("fibpos1Sec2", (fV0AR6-fV0AOctH2+fV0AR5)*cos225/2. - 3.8, (fV0AR6-fV0AOctH2+fV0AR5)*sin225/2.-1.8, 0, rot);
    fibpos1Sec2->RegisterYourself();
    TGeoCompositeShape *sV0AFib1HoleSec2 = new TGeoCompositeShape("sV0AFib1HoleSec2","sV0AFib1Sec2:fibpos1Sec2-sV0AFicR5Sec2"); 
    TGeoVolume *v0AFib1HoleSec2 = new TGeoVolume("V0AFib1HoleSec2",sV0AFib1HoleSec2,medV0AFib);
    v0AFib1HoleSec2->SetLineColor(kV0AColorFib);
    new TGeoArb8("sV0AFib2Sec2", 11.5, v0APts);
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateY(180);
    rot->RotateZ(-90.+22.5);
    TGeoCombiTrans *fibpos2Sec2 = new TGeoCombiTrans("fibpos2Sec2", (fV0AR6-fV0AOctH2+fV0AR5)*cos225/2. - 3.8, (fV0AR6-fV0AOctH2+fV0AR5)*sin225/2.-1.8, 0, rot);
    fibpos2Sec2->RegisterYourself();
    TGeoCompositeShape *sV0AFib2HoleSec2 = new TGeoCompositeShape("sV0AFib2HoleSec2","sV0AFib2Sec2:fibpos2Sec2-sV0AFicR5Sec2");
    TGeoVolume *v0AFib2HoleSec2 = new TGeoVolume("V0AFib2HoleSec2",sV0AFib2HoleSec2,medV0AFib);
    v0AFib2HoleSec2->SetLineColor(kV0AColorFib);
    TGeoVolume *v0AFibSec2 = new TGeoVolumeAssembly("V0AFibSec2");
    v0AFibSec2->AddNode(v0AFib1HoleSec2,1);
    v0AFibSec2->AddNode(v0AFib2HoleSec2,1);
    v0ASec2->AddNode(v0AFibSec2,1); 
    
     /// Plates
    new TGeoTube("sV0ANail1PlaInHoleSec2", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoTube("sV0ANail2PlaInHoleSec2", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoCompositeShape("sV0ANailsPlaInHolesSec2","sV0ANail1PlaInHoleSec2:pos1Sec2+sV0ANail2PlaInHoleSec2:pos2Sec2");
    new TGeoTube("sV0ANail1PlaOuHoleSec2", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoTube("sV0ANail2PlaOuHoleSec2", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoCompositeShape("sV0ANailsPlaOuHolesSec2","sV0ANail1PlaOuHoleSec2:pos1Sec2+sV0ANail2PlaOuHoleSec2:pos2Sec2");
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0;			v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = fV0AR0*sin45;		v0APts[3+8*i] = fV0AR0*sin45;
      v0APts[4+8*i] = fV0AR6 * sin45;	v0APts[5+8*i] = fV0AR6*sin45;
      v0APts[6+8*i] = fV0AR6;		v0APts[7+8*i] = 0.;
    }
    new TGeoArb8("sV0APlaInSec2", (fV0APlaWd-2*fV0APlaAl)/2., v0APts);
    TGeoCompositeShape *sV0APlaInNailsHolesSec2 = new TGeoCompositeShape("sV0APlaInNailsHolesSec2","sV0APlaInSec2-sV0ANailsPlaInHolesSec2-sV0AFicFEEBSec2:posFicFEEBSec2");
    TGeoVolume *v0APlaInNailsHolesSec2 = new TGeoVolume("V0APlaInNailsHolesSec2", sV0APlaInNailsHolesSec2, medV0APlaIn);
    new TGeoArb8("sV0APlaOuSec2", fV0APlaAl/2., v0APts);
    TGeoCompositeShape *sV0APlaOuNailsHolesSec2 = new TGeoCompositeShape("sV0APlaOuNailsHolesSec2","sV0APlaOuSec2-sV0ANailsPlaOuHolesSec2-sV0AFicFEEBSec2:posFicFEEBSec2"); 
    TGeoVolume *v0APlaOuNailsHolesSec2 = new TGeoVolume("V0APlaOuNailsHolesSec2", sV0APlaOuNailsHolesSec2, medV0APlaOu);
    v0APlaInNailsHolesSec2->SetLineColor(kV0AColorPlaIn); v0APlaOuNailsHolesSec2->SetLineColor(kV0AColorPlaOu);
    TGeoVolume *v0APlaSec2 = new TGeoVolumeAssembly("V0APlaSec2");
    v0APlaSec2->AddNode(v0APlaInNailsHolesSec2,1);
    v0APlaSec2->AddNode(v0APlaOuNailsHolesSec2,1,new TGeoTranslation(0,0,(fV0APlaWd-fV0APlaAl)/2.));
    v0APlaSec2->AddNode(v0APlaOuNailsHolesSec2,2,new TGeoTranslation(0,0,-(fV0APlaWd-fV0APlaAl)/2.));
    v0ASec2->AddNode(v0APlaSec2,1,new TGeoTranslation(0,0,(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    v0ASec2->AddNode(v0APlaSec2,2,new TGeoTranslation(0,0,-(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    
     /// Non-sensitive scintilator
    new TGeoTubeSeg("sV0AR5S2Sec2", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShapeSec2, 0, 45);
    TGeoCompositeShape *sV0AR5Sec2 = new TGeoCompositeShape("V0AR5Sec2","sV0AR5S2Sec2 - sV0AChaSec2");
    TGeoVolume *v0AR5Sec2 = new TGeoVolume("V0AR5Sec2",sV0AR5Sec2,medV0ASci);
    v0AR5Sec2->SetLineColor(kV0AColorSci);
    v0ASciSec2->AddNode(v0AR5Sec2,1);
    v0ASec2->AddNode(v0ASciSec2,1); 

    /// PMBox
    TGeoVolume* v0APMSec2 = new TGeoVolumeAssembly("V0APMSec2");
    new TGeoBBox("sV0APMB1Sec2", fV0APMBWd/2., fV0APMBHt/2., fV0APMBTh/2.);
    new TGeoBBox("sV0APMB2Sec2", fV0APMBWd/2.-fV0APMBWdW, fV0APMBHt/2.-fV0APMBHtW, fV0APMBTh/2.-fV0APMBThW);
    TGeoCompositeShape *sV0APMBSec2 = new TGeoCompositeShape("sV0APMBSec2","sV0APMB1Sec2-sV0APMB2Sec2");
    TGeoVolume *v0APMBSec2 = new TGeoVolume("V0APMBSec2",sV0APMBSec2, medV0APMAlum);
    v0APMBSec2->SetLineColor(kV0AColorPMA);
    v0APMSec2->AddNode(v0APMBSec2,1);

    /// PMTubes
    TGeoTube *sV0APMT1Sec2 = new TGeoTube("sV0APMT1Sec2", fV0APMTR1, fV0APMTR2, fV0APMTH/2.);
    TGeoVolume *v0APMT1Sec2 = new TGeoVolume("V0APMT1Sec2", sV0APMT1Sec2, medV0APMGlass);
    TGeoTube *sV0APMT2Sec2 = new TGeoTube("sV0APMT2Sec2", fV0APMTR3, fV0APMTR4, fV0APMTH/2.);
    TGeoVolume *v0APMT2Sec2 = new TGeoVolume("V0APMT2Sec2", sV0APMT2Sec2, medV0APMAlum);
    TGeoVolume *v0APMTSec2 = new TGeoVolumeAssembly("V0APMTSec2");
    TGeoTube *sV0APMTTSec2 = new TGeoTube("sV0APMTTSec2", 0., fV0APMTR4, fV0APMTB/2.);
    TGeoVolume *v0APMTTSec2 = new TGeoVolume("V0APMTTSec2", sV0APMTTSec2, medV0APMAlum);
    v0APMT1Sec2->SetLineColor(kV0AColorPMG);
    v0APMT2Sec2->SetLineColor(kV0AColorPMA);
    v0APMTTSec2->SetLineColor(kV0AColorPMA);
    rot = new TGeoRotation("rot", 90, 0, 180, 0, 90, 90);
    v0APMTSec2->AddNode(v0APMT1Sec2,1,rot);
    v0APMTSec2->AddNode(v0APMT2Sec2,1,rot);
    v0APMTSec2->AddNode(v0APMTTSec2,1,new TGeoCombiTrans(0,-(fV0APMTH+fV0APMTB)/2.,0,rot));
    double autoShiftSec2 = (fV0APMBWd-2*fV0APMBWdW)/4.;
    v0APMSec2->AddNode(v0APMTSec2, 1, new TGeoTranslation(-1.5*autoShiftSec2, 0, 0));
    v0APMSec2->AddNode(v0APMTSec2, 2, new TGeoTranslation(-0.5*autoShiftSec2, 0, 0));
    v0APMSec2->AddNode(v0APMTSec2, 3, new TGeoTranslation(+0.5*autoShiftSec2, 0, 0));
    v0APMSec2->AddNode(v0APMTSec2, 4, new TGeoTranslation(+1.5*autoShiftSec2, 0, 0));

    // PM
    rot = new TGeoRotation("rot");
    rot->RotateX(90-fV0APMBAng);
    rot->RotateZ(-90.+22.5);
    double cosAngPMBSec2 = TMath::Cos(fV0APMBAng*TMath::DegToRad());
    double sinAngPMBSec2 = TMath::Sin(fV0APMBAng*TMath::DegToRad());
    double shiftZSec2 = fV0APMBHt/2. * cosAngPMBSec2
      -   ( fV0ASciWd + 2 * fV0AOctWd + 2 * fV0APlaWd )/2.   -   fV0APMBTh/2. * sinAngPMBSec2;
    double shiftRSec2 = fV0AR6  +  fV0AOctH2 + fV0APlaAl;
    v0ASec2->AddNode(v0APMSec2,1, new TGeoCombiTrans( shiftRSec2*cos225+1.07, shiftRSec2*sin225, shiftZSec2, rot));
    
    // Aluminium nails 
    TGeoTube *sV0ANail1Sec2 = new TGeoTube("sV0ANail1Sec2", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail1Sec2 = new TGeoVolume("V0ANail1Sec2", sV0ANail1Sec2, medV0APMAlum);
    v0ANail1Sec2->SetLineColor(kV0AColorPMA);// this is the color for aluminium
    v0ASec2->AddNode(v0ANail1Sec2,1,new TGeoTranslation(42.9, 0.51, 0.0));
    TGeoTube *sV0ANail2Sec2 = new TGeoTube("sV0ANail2Sec2", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail2Sec2 = new TGeoVolume("V0ANail2Sec2", sV0ANail2Sec2, medV0APMAlum);
    v0ANail2Sec2->SetLineColor(kV0AColorPMA);
    v0ASec2->AddNode(v0ANail2Sec2,1,new TGeoTranslation(30.73,29.98,0.0)); 
        
    /// Adding sector to v0LE volume
    for(int i=1; i<2; i++) {
       TGeoRotation *rotation = new TGeoRotation("rotation", 90., i*45., 90., 90.+i*45., 0., 0.);
       v0LE->AddNode(v0ASec2,i+1,rotation);  
    }

    //FEEBox
    TGeoVolume* v0AFEE2 = new TGeoVolumeAssembly("V0AFEE2");
    v0AFEE2->AddNode(v0AFEEB,1);

    //Mother and daughter boards
    v0AFEE2->AddNode(v0AFEEDaughter, 1, new TGeoTranslation(9.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE2->AddNode(v0AFEEDaughter, 2, new TGeoTranslation(6.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE2->AddNode(v0AFEEDaughter, 3, new TGeoTranslation(3.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE2->AddNode(v0AFEEDaughter, 4, new TGeoTranslation(0.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE2->AddNode(v0AFEEDaughter, 5, new TGeoTranslation(-3.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE2->AddNode(v0AFEEDaughter, 6, new TGeoTranslation(-6.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE2->AddNode(v0AFEEDaughter, 7, new TGeoTranslation(-9.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE2->AddNode(v0AFEEDaughter, 8, new TGeoTranslation(-12.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));   
    v0AFEE2->AddNode(v0AFEEMother, 1, new TGeoTranslation(0.0, 0.0, -fV0AFEEBTh/2.+fV0APMBThW+fV0APMBHtW));
    v0AFEE2->AddNode(v0AFEEHalfMother, 1, new TGeoTranslation(0.0, -(fV0AFEEBTh+fV0APMTB)/2., -2.0*spacing));

    //FEE
    rot = new TGeoRotation("rot");
    rot->RotateX(90);
    rot->RotateZ(+90.0-45.0+90);
    v0LE->AddNode(v0AFEE2,1, new TGeoCombiTrans( (aFEEshiftR2Sec1*cos225 + 2.0)*cos45, (aFEEshiftR2Sec1*cos225 + 2.0)*sin45, 7.5, rot));
    rot = new TGeoRotation("rot");
    rot->RotateZ(-90.0+45.0+90.0);
    v0LE->AddNode(v0AFEEOct2,2, new TGeoCombiTrans(((aFEEshiftR2Sec1*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0)*cos45,((aFEEshiftR2Sec1*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0)*sin45,0.0,rot));

   
    //Upper supports
    for (int i=0;i<2;i++){
    v0APts[0+8*i] = 0.0;	    v0APts[1+8*i] = 45.5;  
    v0APts[2+8*i] = 0.0;            v0APts[3+8*i] = 70.4;   
    v0APts[4+8*i] = 4.0;	    v0APts[5+8*i] = 68.9;
    v0APts[6+8*i] = 4.0;	    v0APts[7+8*i] = 45.5;  
    }
    TGeoArb8 *sV0ASuppur = new TGeoArb8("sV0ASuppur", 1.65, v0APts);    
    TGeoVolume *v0ASuppur = new TGeoVolume("V0ASuppur", sV0ASuppur, medV0ASup);
    v0ASuppur->SetLineColor(kV0AColorOct);
    v0LE->AddNode(v0ASuppur,1);
    for (int i=0;i<2;i++){
    v0APts[0+8*i] = -0.0;	    v0APts[1+8*i] = 70.4;
    v0APts[2+8*i] = -0.0;           v0APts[3+8*i] = 45.5;
    v0APts[4+8*i] = -4.0;	    v0APts[5+8*i] = 45.5;
    v0APts[6+8*i] = -4.0;	    v0APts[7+8*i] = 68.9;
    }
    TGeoArb8 *sV0ASuppul = new TGeoArb8("sV0ASuppul", 1.65, v0APts);    
    TGeoVolume *v0ASuppul = new TGeoVolume("V0ASuppul", sV0ASuppul, medV0ASup);
    v0ASuppul->SetLineColor(kV0AColorOct);
    v0LE->AddNode(v0ASuppul,1);

    /// Definition sector 3
    TGeoVolume *v0ASec3 = new TGeoVolumeAssembly("V0ASec3");
        
    /// For boolean sustraction
    double preShapeSec3 = 0.2;
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0-fV0AFraWd/2.-preShapeSec3;  v0APts[1+8*i] = -preShapeSec3;
      v0APts[2+8*i] = fV0AR0-fV0AFraWd/2.-preShapeSec3;  v0APts[3+8*i] = fV0AFraWd/2.;
      v0APts[4+8*i] = fV0AR4+fV0AFraWd/2.+preShapeSec3;  v0APts[5+8*i] = fV0AFraWd/2.;
      v0APts[6+8*i] = fV0AR4+fV0AFraWd/2.+preShapeSec3;  v0APts[7+8*i] = -preShapeSec3;
    }
    new TGeoArb8("sV0ACha1Sec3",fV0ASciWd/1.5,v0APts);
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0*sin45-preShapeSec3;
      v0APts[1+8*i] = (fV0AR0-fV0AFraWd)*sin45-preShapeSec3;
      v0APts[2+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45-preShapeSec3;
      v0APts[3+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45;
      v0APts[4+8*i] = (fV0AR4+fV0AFraWd/2.)*sin45+preShapeSec3;
      v0APts[5+8*i] = (fV0AR4+fV0AFraWd/2.)*sin45+2.*preShapeSec3;
      v0APts[6+8*i] = (fV0AR4+fV0AFraWd)*sin45+preShapeSec3;
      v0APts[7+8*i] = fV0AR4*sin45+preShapeSec3;
    }
    new TGeoArb8("sV0ACha2Sec3", fV0ASciWd/2.+2.*preShapeSec3, v0APts);
    new TGeoCompositeShape("sV0ACha12Sec3","sV0ACha1Sec3+sV0ACha2Sec3");
    new TGeoTube("sV0ANail1SciHoleSec3", 0.0, 0.4, 1.65);
    TGeoTranslation *pos1Sec3 = new TGeoTranslation("pos1Sec3", 42.9, 0.51, 0.0);
    pos1Sec3->RegisterYourself();
    new TGeoTube("sV0ANail2SciHoleSec3", 0.0, 0.4, 1.65);
    TGeoTranslation *pos2Sec3 = new TGeoTranslation("pos2Sec3", 30.73,29.98,0.0);
    pos2Sec3->RegisterYourself();
    new TGeoCompositeShape("sV0ANailsSciHolesSec3","sV0ANail1SciHoleSec3:pos1Sec3+sV0ANail2SciHoleSec3:pos2Sec3");
    new TGeoCompositeShape("sV0AChaSec3","sV0ACha12Sec3+sV0ANailsSciHolesSec3");
    new TGeoTubeSeg("sV0AFicR5Sec3", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShapeSec3, 0, 45);
    new TGeoBBox("sV0AFicFEEBSec3", fV0AFEEBWd/2., fV0AFEEBHt/2., fV0AFEEBTh/2.);
    double aFEEshiftR2Sec3 = fV0AR6  +  fV0AOctH2 + fV0APlaAl;
    rot = new TGeoRotation("rot");
    rot->RotateX(90);
    rot->RotateZ(-90.0+45.0);
    TGeoCombiTrans *posFicFEEBSec3 = new TGeoCombiTrans("posFicFEEBSec3", (aFEEshiftR2Sec3*cos225 + 2.0 )*cos45, (aFEEshiftR2Sec3*cos225 + 2.0 )*sin45, 7.5, rot);
    posFicFEEBSec3->RegisterYourself();
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = fV0AOctH2/2.;	       v0APts[1+8*i] = fV0AFEEBHt/2. + 2.5;
    v0APts[2+8*i] = fV0AOctH2/2.;              v0APts[3+8*i] = -fV0AFEEBHt/2.- 2.5;
    v0APts[4+8*i] = -fV0AOctH2/2.;	       v0APts[5+8*i] = -fV0AFEEBHt/2.- 2.5;
    v0APts[6+8*i] = -fV0AOctH2/2.;	       v0APts[7+8*i] = fV0AFEEBHt/2.+ 2.5;
    }
    new TGeoArb8("sV0AFicOct2Sec3", (fV0ASciWd+2*fV0AOctWd)/2., v0APts); 
    rot = new TGeoRotation("rot");
    rot->RotateZ(-90.0+45.0+90.0);
    TGeoCombiTrans *posFicOct2UpSec3 = new TGeoCombiTrans("posFicOct2UpSec3",((aFEEshiftR2Sec3*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0)*cos45,((aFEEshiftR2Sec3*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0)*sin45,0.0,rot);
    posFicOct2UpSec3->RegisterYourself();

    /// Frame
    TGeoVolume *v0AFraSec3 = new TGeoVolumeAssembly("V0AFraSec3");
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[3+8*i] = fV0AFraWd/2.;
      v0APts[4+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[5+8*i] = fV0AFraWd/2.;
      v0APts[6+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[7+8*i] = 0.;
    }
    TGeoArb8 *sV0AFraB1Sec3 = new TGeoArb8("sV0AFraB1Sec3",fV0ASciWd/2.,v0APts);
    TGeoVolume *v0AFraB1Sec3 = new TGeoVolume("V0AFraB1Sec3",sV0AFraB1Sec3,medV0AFra);
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
    TGeoArb8 *sV0AFraB2Sec3 = new TGeoArb8("sV0AFraB2Sec3", fV0ASciWd/2., v0APts);
    TGeoVolume *v0AFraB2Sec3 = new TGeoVolume("V0AFraB2Sec3",sV0AFraB2Sec3,medV0AFra);
    v0AFraB1Sec3->SetLineColor(kV0AColorFra); v0AFraB2Sec3->SetLineColor(kV0AColorFra);
    v0AFraSec3->AddNode(v0AFraB1Sec3,1);
    v0AFraSec3->AddNode(v0AFraB2Sec3,1);  
    new TGeoTubeSeg( "sV0AFraR1bSec3", fV0AR0-fV0AFraWd/2.,
		     fV0AR0+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR2bSec3", fV0AR1-fV0AFraWd/2.,
		     fV0AR1+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR3bSec3", fV0AR2-fV0AFraWd/2.,
		     fV0AR2+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR4bSec3", fV0AR3-fV0AFraWd/2.,
		     fV0AR3+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR5bSec3", fV0AR4-fV0AFraWd/2.,
		     fV0AR4+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    TGeoCompositeShape *sV0AFraR1Sec3 = new TGeoCompositeShape("sV0AFraR1Sec3","sV0AFraR1bSec3-sV0AChaSec3");
    TGeoCompositeShape *sV0AFraR2Sec3 = new TGeoCompositeShape("sV0AFraR2Sec3","sV0AFraR2bSec3-sV0AChaSec3");
    TGeoCompositeShape *sV0AFraR3Sec3 = new TGeoCompositeShape("sV0AFraR3Sec3","sV0AFraR3bSec3-sV0AChaSec3");
    TGeoCompositeShape *sV0AFraR4Sec3 = new TGeoCompositeShape("sV0AFraR4Sec3","sV0AFraR4bSec3-sV0AChaSec3");
    TGeoCompositeShape *sV0AFraR5Sec3 = new TGeoCompositeShape("sV0AFraR5Sec3","sV0AFraR5bSec3-sV0AChaSec3");
    TGeoVolume *v0AFraR1Sec3 = new TGeoVolume("V0AFraR1Sec3",sV0AFraR1Sec3,medV0AFra);
    TGeoVolume *v0AFraR2Sec3 = new TGeoVolume("V0AFraR2Sec3",sV0AFraR2Sec3,medV0AFra);
    TGeoVolume *v0AFraR3Sec3 = new TGeoVolume("V0AFraR3Sec3",sV0AFraR3Sec3,medV0AFra);
    TGeoVolume *v0AFraR4Sec3 = new TGeoVolume("V0AFraR4Sec3",sV0AFraR4Sec3,medV0AFra);
    TGeoVolume *v0AFraR5Sec3 = new TGeoVolume("V0AFraR5Sec3",sV0AFraR5Sec3,medV0AFra);
    v0AFraR1Sec3->SetLineColor(kV0AColorFra); v0AFraR2Sec3->SetLineColor(kV0AColorFra);
    v0AFraR3Sec3->SetLineColor(kV0AColorFra); v0AFraR4Sec3->SetLineColor(kV0AColorFra);
    v0AFraR5Sec3->SetLineColor(kV0AColorFra);
    v0AFraSec3->AddNode(v0AFraR1Sec3,1); 
    v0AFraSec3->AddNode(v0AFraR2Sec3,1);
    v0AFraSec3->AddNode(v0AFraR3Sec3,1); 
    v0AFraSec3->AddNode(v0AFraR4Sec3,1);
    v0AFraSec3->AddNode(v0AFraR5Sec3,1);
    v0ASec3->AddNode(v0AFraSec3,1);
    
    /// Sensitive scintilator
    TGeoVolume *v0ASciSec3 = new TGeoVolumeAssembly("V0ASciSec3");
    new TGeoTubeSeg( "sV0AR1bSec3", fV0AR0+fV0AFraWd/2.,
		     fV0AR1-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR2bSec3", fV0AR1+fV0AFraWd/2.,
		     fV0AR2-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR3bSec3", fV0AR2+fV0AFraWd/2.,
		     fV0AR3-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR4bSec3", fV0AR3+fV0AFraWd/2.,
		     fV0AR4-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    TGeoCompositeShape *sV0AR1Sec3 = new TGeoCompositeShape("sV0AR1Sec3","sV0AR1bSec3-sV0AChaSec3");
    TGeoCompositeShape *sV0AR2Sec3 = new TGeoCompositeShape("sV0AR2Sec3","sV0AR2bSec3-sV0AChaSec3");
    TGeoCompositeShape *sV0AR3Sec3 = new TGeoCompositeShape("sV0AR3Sec3","sV0AR3bSec3-sV0AChaSec3");
    TGeoCompositeShape *sV0AR4Sec3 = new TGeoCompositeShape("sV0AR4Sec3","sV0AR4bSec3-sV0AChaSec3");
    TGeoVolume *v0L1Sec3 = new TGeoVolume("V0L1Sec3",sV0AR1Sec3,medV0ASci);
    TGeoVolume *v0L2Sec3 = new TGeoVolume("V0L2Sec3",sV0AR2Sec3,medV0ASci);
    TGeoVolume *v0L3Sec3 = new TGeoVolume("V0L3Sec3",sV0AR3Sec3,medV0ASci);
    TGeoVolume *v0L4Sec3 = new TGeoVolume("V0L4Sec3",sV0AR4Sec3,medV0ASci);
    v0L1Sec3->SetLineColor(kV0AColorSci); v0L2Sec3->SetLineColor(kV0AColorSci);
    v0L3Sec3->SetLineColor(kV0AColorSci); v0L4Sec3->SetLineColor(kV0AColorSci);
    v0ASec3->AddNode(v0L1Sec3,1);
    v0ASec3->AddNode(v0L2Sec3,1);
    v0ASec3->AddNode(v0L3Sec3,1);
    v0ASec3->AddNode(v0L4Sec3,1);      
    
    /// Segment of octagon 
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] =  fV0AR6-fV0AOctH2;           v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = (fV0AR6-fV0AOctH2)*sin45;  v0APts[3+8*i] = (fV0AR6-fV0AOctH2)*sin45;
      v0APts[4+8*i] = fV0AR6*sin45;		 v0APts[5+8*i] = fV0AR6*sin45;
      v0APts[6+8*i] = fV0AR6;			 v0APts[7+8*i] = 0.;
    }
    new TGeoArb8("sV0AOct2Sec3", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoCompositeShape *sV0AOct2FEEBSec3 = new TGeoCompositeShape("sV0AOct2FEEBSec3","sV0AOct2Sec3-sV0AFicOct2Sec3:posFicOct2UpSec3-sV0AFicFEEBSec3:posFicFEEBSec3");
    TGeoVolume *v0AOct2Sec3 = new TGeoVolume("V0AOct2Sec3", sV0AOct2FEEBSec3,medV0ASup);
    v0AOct2Sec3->SetLineColor(kV0AColorOct);
    TGeoVolume *v0ASupSec3 = new TGeoVolumeAssembly("V0ASupSec3");
    v0ASupSec3->AddNode(v0AOct2Sec3,1);
    v0ASec3->AddNode(v0ASupSec3,1);

    //Bunch of fibers
    v0APts[ 0] = v0APts[ 2] = -13.0;
    v0APts[ 1] = v0APts[ 7] = (fV0ASciWd+fV0AOctWd)/2.-0.01;
    v0APts[ 3] = v0APts[ 5] = (fV0ASciWd+fV0AOctWd)/2.+0.01;
    v0APts[ 4] = v0APts[ 6] = +13.0;
    v0APts[ 8] = v0APts[10] = -10.0;
    v0APts[ 9] = v0APts[15] = 0.;
    v0APts[11] = v0APts[13] = 0.25;
    v0APts[12] = v0APts[14] = +10.0;
    new TGeoArb8("sV0AFib1Sec3", 11.5, v0APts);
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateZ(-90.+22.5);
    TGeoCombiTrans *fibpos1Sec3 = new TGeoCombiTrans("fibpos1Sec3", (fV0AR6-fV0AOctH2+fV0AR5)*cos225/2. - 3.8, (fV0AR6-fV0AOctH2+fV0AR5)*sin225/2.-1.8, 0, rot);
    fibpos1Sec3->RegisterYourself();
    TGeoCompositeShape *sV0AFib1HoleSec3 = new TGeoCompositeShape("sV0AFib1HoleSec3","sV0AFib1Sec3:fibpos1Sec3-sV0AFicR5Sec3"); 
    TGeoVolume *v0AFib1HoleSec3 = new TGeoVolume("V0AFib1HoleSec3",sV0AFib1HoleSec3,medV0AFib);
    v0AFib1HoleSec3->SetLineColor(kV0AColorFib);
    new TGeoArb8("sV0AFib2Sec3", 11.5, v0APts);
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateY(180);
    rot->RotateZ(-90.+22.5);
    TGeoCombiTrans *fibpos2Sec3 = new TGeoCombiTrans("fibpos2Sec3", (fV0AR6-fV0AOctH2+fV0AR5)*cos225/2. - 3.8, (fV0AR6-fV0AOctH2+fV0AR5)*sin225/2.-1.8, 0, rot);
    fibpos2Sec3->RegisterYourself();
    TGeoCompositeShape *sV0AFib2HoleSec3 = new TGeoCompositeShape("sV0AFib2HoleSec3","sV0AFib2Sec3:fibpos2Sec3-sV0AFicR5Sec3");
    TGeoVolume *v0AFib2HoleSec3 = new TGeoVolume("V0AFib2HoleSec3",sV0AFib2HoleSec3,medV0AFib);
    v0AFib2HoleSec3->SetLineColor(kV0AColorFib);
    TGeoVolume *v0AFibSec3 = new TGeoVolumeAssembly("V0AFibSec3");
    v0AFibSec3->AddNode(v0AFib1HoleSec3,1);
    v0AFibSec3->AddNode(v0AFib2HoleSec3,1);
    v0ASec3->AddNode(v0AFibSec3,1); 
    
     /// Plates
    new TGeoTube("sV0ANail1PlaInHoleSec3", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoTube("sV0ANail2PlaInHoleSec3", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoCompositeShape("sV0ANailsPlaInHolesSec3","sV0ANail1PlaInHoleSec3:pos1Sec3+sV0ANail2PlaInHoleSec3:pos2Sec3");
    new TGeoTube("sV0ANail1PlaOuHoleSec3", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoTube("sV0ANail2PlaOuHoleSec3", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoCompositeShape("sV0ANailsPlaOuHolesSec3","sV0ANail1PlaOuHoleSec3:pos1Sec3+sV0ANail2PlaOuHoleSec3:pos2Sec3");
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0;			v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = fV0AR0*sin45;		v0APts[3+8*i] = fV0AR0*sin45;
      v0APts[4+8*i] = fV0AR6 * sin45;	v0APts[5+8*i] = fV0AR6*sin45;
      v0APts[6+8*i] = fV0AR6;		v0APts[7+8*i] = 0.;
    }
    new TGeoArb8("sV0APlaInSec3", (fV0APlaWd-2*fV0APlaAl)/2., v0APts);
    TGeoCompositeShape *sV0APlaInNailsHolesSec3 = new TGeoCompositeShape("sV0APlaInNailsHolesSec3","sV0APlaInSec3-sV0ANailsPlaInHolesSec3-sV0AFicFEEBSec3:posFicFEEBSec3");
    TGeoVolume *v0APlaInNailsHolesSec3 = new TGeoVolume("V0APlaInNailsHolesSec3", sV0APlaInNailsHolesSec3, medV0APlaIn);
    new TGeoArb8("sV0APlaOuSec3", fV0APlaAl/2., v0APts);
    TGeoCompositeShape *sV0APlaOuNailsHolesSec3 = new TGeoCompositeShape("sV0APlaOuNailsHolesSec3","sV0APlaOuSec3-sV0ANailsPlaOuHolesSec3-sV0AFicFEEBSec3:posFicFEEBSec3"); 
    TGeoVolume *v0APlaOuNailsHolesSec3 = new TGeoVolume("V0APlaOuNailsHolesSec3", sV0APlaOuNailsHolesSec3, medV0APlaOu);
    v0APlaInNailsHolesSec3->SetLineColor(kV0AColorPlaIn); v0APlaOuNailsHolesSec3->SetLineColor(kV0AColorPlaOu);
    TGeoVolume *v0APlaSec3 = new TGeoVolumeAssembly("V0APlaSec3");
    v0APlaSec3->AddNode(v0APlaInNailsHolesSec3,1);
    v0APlaSec3->AddNode(v0APlaOuNailsHolesSec3,1,new TGeoTranslation(0,0,(fV0APlaWd-fV0APlaAl)/2.));
    v0APlaSec3->AddNode(v0APlaOuNailsHolesSec3,2,new TGeoTranslation(0,0,-(fV0APlaWd-fV0APlaAl)/2.));
    v0ASec3->AddNode(v0APlaSec3,1,new TGeoTranslation(0,0,(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    v0ASec3->AddNode(v0APlaSec3,2,new TGeoTranslation(0,0,-(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    
     /// Non-sensitive scintilator
    new TGeoTubeSeg("sV0AR5S2Sec3", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShapeSec3, 0, 45);
    TGeoCompositeShape *sV0AR5Sec3 = new TGeoCompositeShape("V0AR5Sec3","sV0AR5S2Sec3 - sV0AChaSec3");
    TGeoVolume *v0AR5Sec3 = new TGeoVolume("V0AR5Sec3",sV0AR5Sec3,medV0ASci);
    v0AR5Sec3->SetLineColor(kV0AColorSci);
    v0ASciSec3->AddNode(v0AR5Sec3,1);
    v0ASec3->AddNode(v0ASciSec3,1); 

    /// PMBox
    TGeoVolume* v0APMSec3 = new TGeoVolumeAssembly("V0APMSec3");
    new TGeoBBox("sV0APMB1Sec3", fV0APMBWd/2., fV0APMBHt/2., fV0APMBTh/2.);
    new TGeoBBox("sV0APMB2Sec3", fV0APMBWd/2.-fV0APMBWdW, fV0APMBHt/2.-fV0APMBHtW, fV0APMBTh/2.-fV0APMBThW);
    TGeoCompositeShape *sV0APMBSec3 = new TGeoCompositeShape("sV0APMBSec3","sV0APMB1Sec3-sV0APMB2Sec3");
    TGeoVolume *v0APMBSec3 = new TGeoVolume("V0APMBSec3",sV0APMBSec3, medV0APMAlum);
    v0APMBSec3->SetLineColor(kV0AColorPMA);
    v0APMSec3->AddNode(v0APMBSec3,1);

    /// PMTubes
    TGeoTube *sV0APMT1Sec3 = new TGeoTube("sV0APMT1Sec3", fV0APMTR1, fV0APMTR2, fV0APMTH/2.);
    TGeoVolume *v0APMT1Sec3 = new TGeoVolume("V0APMT1Sec3", sV0APMT1Sec3, medV0APMGlass);
    TGeoTube *sV0APMT2Sec3 = new TGeoTube("sV0APMT2Sec3", fV0APMTR3, fV0APMTR4, fV0APMTH/2.);
    TGeoVolume *v0APMT2Sec3 = new TGeoVolume("V0APMT2Sec3", sV0APMT2Sec3, medV0APMAlum);
    TGeoVolume *v0APMTSec3 = new TGeoVolumeAssembly("V0APMTSec3");
    TGeoTube *sV0APMTTSec3 = new TGeoTube("sV0APMTTSec3", 0., fV0APMTR4, fV0APMTB/2.);
    TGeoVolume *v0APMTTSec3 = new TGeoVolume("V0APMTTSec3", sV0APMTTSec3, medV0APMAlum);
    v0APMT1Sec3->SetLineColor(kV0AColorPMG);
    v0APMT2Sec3->SetLineColor(kV0AColorPMA);
    v0APMTTSec3->SetLineColor(kV0AColorPMA);
    rot = new TGeoRotation("rot", 90, 0, 180, 0, 90, 90);
    v0APMTSec3->AddNode(v0APMT1Sec3,1,rot);
    v0APMTSec3->AddNode(v0APMT2Sec3,1,rot);
    v0APMTSec3->AddNode(v0APMTTSec3,1,new TGeoCombiTrans(0,-(fV0APMTH+fV0APMTB)/2.,0,rot));
    double autoShiftSec3 = (fV0APMBWd-2*fV0APMBWdW)/4.;
    v0APMSec3->AddNode(v0APMTSec3, 1, new TGeoTranslation(-1.5*autoShiftSec3, 0, 0));
    v0APMSec3->AddNode(v0APMTSec3, 2, new TGeoTranslation(-0.5*autoShiftSec3, 0, 0));
    v0APMSec3->AddNode(v0APMTSec3, 3, new TGeoTranslation(+0.5*autoShiftSec3, 0, 0));
    v0APMSec3->AddNode(v0APMTSec3, 4, new TGeoTranslation(+1.5*autoShiftSec3, 0, 0));

    // PM
    rot = new TGeoRotation("rot");
    rot->RotateX(90-fV0APMBAng);
    rot->RotateZ(-90.+22.5);
    double cosAngPMBSec3 = TMath::Cos(fV0APMBAng*TMath::DegToRad());
    double sinAngPMBSec3 = TMath::Sin(fV0APMBAng*TMath::DegToRad());
    double shiftZSec3 = fV0APMBHt/2. * cosAngPMBSec3
      -   ( fV0ASciWd + 2 * fV0AOctWd + 2 * fV0APlaWd )/2.   -   fV0APMBTh/2. * sinAngPMBSec3;
    double shiftRSec3 = fV0AR6  +  fV0AOctH2 + fV0APlaAl;
    v0ASec3->AddNode(v0APMSec3,1, new TGeoCombiTrans( shiftRSec3*cos225+1.07, shiftRSec3*sin225, shiftZSec3, rot));
    
    // Aluminium nails 
    TGeoTube *sV0ANail1Sec3 = new TGeoTube("sV0ANail1Sec3", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail1Sec3 = new TGeoVolume("V0ANail1Sec3", sV0ANail1Sec3, medV0APMAlum);
    v0ANail1Sec3->SetLineColor(kV0AColorPMA);// this is the color for aluminium
    v0ASec3->AddNode(v0ANail1Sec3,1,new TGeoTranslation(42.9, 0.51, 0.0));
    TGeoTube *sV0ANail2Sec3 = new TGeoTube("sV0ANail2Sec3", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail2Sec3 = new TGeoVolume("V0ANail2Sec3", sV0ANail2Sec3, medV0APMAlum);
    v0ANail2Sec3->SetLineColor(kV0AColorPMA);
    v0ASec3->AddNode(v0ANail2Sec3,1,new TGeoTranslation(30.73,29.98,0.0)); 
        
    /// Adding sector to v0LE volume
    for(int i=2; i<3; i++) {
       TGeoRotation *rotation = new TGeoRotation("rotation", 90., i*45., 90., 90.+i*45., 0., 0.);
       v0LE->AddNode(v0ASec3,i+1,rotation);  
    }

    //FEEBox
    TGeoVolume* v0AFEE3 = new TGeoVolumeAssembly("V0AFEE3");
    v0AFEE3->AddNode(v0AFEEB,1);

    //Mother and daughter boards
    v0AFEE3->AddNode(v0AFEEDaughter, 1, new TGeoTranslation(9.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE3->AddNode(v0AFEEDaughter, 2, new TGeoTranslation(6.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE3->AddNode(v0AFEEDaughter, 3, new TGeoTranslation(3.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE3->AddNode(v0AFEEDaughter, 4, new TGeoTranslation(0.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE3->AddNode(v0AFEEDaughter, 5, new TGeoTranslation(-3.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE3->AddNode(v0AFEEDaughter, 6, new TGeoTranslation(-6.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE3->AddNode(v0AFEEDaughter, 7, new TGeoTranslation(-9.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE3->AddNode(v0AFEEDaughter, 8, new TGeoTranslation(-12.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));   
    v0AFEE3->AddNode(v0AFEEMother, 1, new TGeoTranslation(0.0, 0.0, -fV0AFEEBTh/2.+fV0APMBThW+fV0APMBHtW));
    v0AFEE3->AddNode(v0AFEEHalfMother, 1, new TGeoTranslation(0.0, -(fV0AFEEBTh+fV0APMTB)/2., -2.0*spacing));

    //FEE
    rot = new TGeoRotation("rot");
    rot->RotateX(90);
    rot->RotateZ(-90.0+45.0-90);
    v0LE->AddNode(v0AFEE3,1, new TGeoCombiTrans( -(aFEEshiftR2Sec1*cos225 + 2.0)*cos45, (aFEEshiftR2Sec1*cos225 + 2.0)*sin45, 7.5, rot) );
    rot = new TGeoRotation("rot");
    rot->RotateZ(+90.0-45.0-90.0);
    v0LE->AddNode(v0AFEEOct2,3, new TGeoCombiTrans(-((aFEEshiftR2Sec1*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0)*cos45,((aFEEshiftR2Sec1*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0)*sin45,0.0,rot)); 


    /// Definition sector 4
    TGeoVolume *v0ASec4 = new TGeoVolumeAssembly("V0ASec4");
        
    /// For boolean sustraction
    double preShapeSec4 = 0.2;
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0-fV0AFraWd/2.-preShapeSec4;  v0APts[1+8*i] = -preShapeSec4;
      v0APts[2+8*i] = fV0AR0-fV0AFraWd/2.-preShapeSec4;  v0APts[3+8*i] = fV0AFraWd/2.;
      v0APts[4+8*i] = fV0AR4+fV0AFraWd/2.+preShapeSec4;  v0APts[5+8*i] = fV0AFraWd/2.;
      v0APts[6+8*i] = fV0AR4+fV0AFraWd/2.+preShapeSec4;  v0APts[7+8*i] = -preShapeSec4;
    }
    new TGeoArb8("sV0ACha1Sec4",fV0ASciWd/1.5,v0APts);
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0*sin45-preShapeSec4;
      v0APts[1+8*i] = (fV0AR0-fV0AFraWd)*sin45-preShapeSec4;
      v0APts[2+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45-preShapeSec4;
      v0APts[3+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45;
      v0APts[4+8*i] = (fV0AR4+fV0AFraWd/2.)*sin45+preShapeSec4;
      v0APts[5+8*i] = (fV0AR4+fV0AFraWd/2.)*sin45+2.*preShapeSec4;
      v0APts[6+8*i] = (fV0AR4+fV0AFraWd)*sin45+preShapeSec4;
      v0APts[7+8*i] = fV0AR4*sin45+preShapeSec4;
    }
    new TGeoArb8("sV0ACha2Sec4", fV0ASciWd/2.+2.*preShapeSec4, v0APts);
    new TGeoCompositeShape("sV0ACha12Sec4","sV0ACha1Sec4+sV0ACha2Sec4");
    new TGeoTube("sV0ANail1SciHoleSec4", 0.0, 0.4, 1.65);
    TGeoTranslation *pos1Sec4 = new TGeoTranslation("pos1Sec4", 42.9, 0.51, 0.0);
    pos1Sec4->RegisterYourself();
    new TGeoTube("sV0ANail2SciHoleSec4", 0.0, 0.4, 1.65);
    TGeoTranslation *pos2Sec4 = new TGeoTranslation("pos2Sec4", 30.73,29.98,0.0);
    pos2Sec4->RegisterYourself();
    new TGeoCompositeShape("sV0ANailsSciHolesSec4","sV0ANail1SciHoleSec4:pos1Sec4+sV0ANail2SciHoleSec4:pos2Sec4");
    new TGeoCompositeShape("sV0AChaSec4","sV0ACha12Sec4+sV0ANailsSciHolesSec4");
    new TGeoTubeSeg("sV0AFicR5Sec4", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShapeSec4, 0, 45);
    new TGeoBBox("sV0AFicFEEBSec4", fV0AFEEBWd/2., fV0AFEEBHt/2., fV0AFEEBTh/2.);
    rot = new TGeoRotation("rot");
    rot->RotateX(90);
    rot->RotateZ(-90.0);
    double aFEEshiftR2Sec4 = fV0AR6  +  fV0AOctH2 + fV0APlaAl;
    TGeoCombiTrans *posFicFEEBSec4 = new TGeoCombiTrans("posFicFEEBSec4", aFEEshiftR2Sec4*cos225 + 2.0, 0, 7.5, rot);
    posFicFEEBSec4->RegisterYourself();
    rot = new TGeoRotation("rot");
    rot->RotateX(90);
    rot->RotateZ(-90.0+45.0);
    TGeoCombiTrans *posFicFEEBUpSec4 = new TGeoCombiTrans("posFicFEEBUpSec4", (aFEEshiftR2Sec4*cos225 + 2.0 )*cos45, (aFEEshiftR2Sec4*cos225 + 2.0 )*sin45, 7.5, rot);
    posFicFEEBUpSec4->RegisterYourself();
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = fV0AOctH2/2.;	       v0APts[1+8*i] = fV0AFEEBHt/2. + 2.5;
    v0APts[2+8*i] = fV0AOctH2/2.;              v0APts[3+8*i] = -fV0AFEEBHt/2.- 2.5;
    v0APts[4+8*i] = -fV0AOctH2/2.;	       v0APts[5+8*i] = -fV0AFEEBHt/2.- 2.5;
    v0APts[6+8*i] = -fV0AOctH2/2.;	       v0APts[7+8*i] = fV0AFEEBHt/2.+ 2.5;
    }
    new TGeoArb8("sV0AFicOct2Sec4", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoTranslation *posFicOct2Sec4 = new TGeoTranslation("posFicOct2Sec4",(aFEEshiftR2Sec4*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0,0.0,0.0);
    posFicOct2Sec4->RegisterYourself();  
    rot = new TGeoRotation("rot");
    rot->RotateZ(-90.0+45.0+90.0);
    TGeoCombiTrans *posFicOct2UpSec4 = new TGeoCombiTrans("posFicOct2UpSec4",((aFEEshiftR2Sec4*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0)*cos45,((aFEEshiftR2Sec4*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0)*sin45,0.0,rot);
    posFicOct2UpSec4->RegisterYourself(); 

    /// Frame
    TGeoVolume *v0AFraSec4 = new TGeoVolumeAssembly("V0AFraSec4");
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[3+8*i] = fV0AFraWd/2.;
      v0APts[4+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[5+8*i] = fV0AFraWd/2.;
      v0APts[6+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[7+8*i] = 0.;
    }
    TGeoArb8 *sV0AFraB1Sec4 = new TGeoArb8("sV0AFraB1Sec4",fV0ASciWd/2.,v0APts);
    TGeoVolume *v0AFraB1Sec4 = new TGeoVolume("V0AFraB1Sec4",sV0AFraB1Sec4,medV0AFra);
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
    TGeoArb8 *sV0AFraB2Sec4 = new TGeoArb8("sV0AFraB2Sec4", fV0ASciWd/2., v0APts);
    TGeoVolume *v0AFraB2Sec4 = new TGeoVolume("V0AFraB2Sec4",sV0AFraB2Sec4,medV0AFra);
    v0AFraB1Sec4->SetLineColor(kV0AColorFra); v0AFraB2Sec4->SetLineColor(kV0AColorFra);
    v0AFraSec4->AddNode(v0AFraB1Sec4,1);
    v0AFraSec4->AddNode(v0AFraB2Sec4,1);  
    new TGeoTubeSeg( "sV0AFraR1bSec4", fV0AR0-fV0AFraWd/2.,
		     fV0AR0+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR2bSec4", fV0AR1-fV0AFraWd/2.,
		     fV0AR1+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR3bSec4", fV0AR2-fV0AFraWd/2.,
		     fV0AR2+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR4bSec4", fV0AR3-fV0AFraWd/2.,
		     fV0AR3+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AFraR5bSec4", fV0AR4-fV0AFraWd/2.,
		     fV0AR4+fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    TGeoCompositeShape *sV0AFraR1Sec4 = new TGeoCompositeShape("sV0AFraR1Sec4","sV0AFraR1bSec4-sV0AChaSec4");
    TGeoCompositeShape *sV0AFraR2Sec4 = new TGeoCompositeShape("sV0AFraR2Sec4","sV0AFraR2bSec4-sV0AChaSec4");
    TGeoCompositeShape *sV0AFraR3Sec4 = new TGeoCompositeShape("sV0AFraR3Sec4","sV0AFraR3bSec4-sV0AChaSec4");
    TGeoCompositeShape *sV0AFraR4Sec4 = new TGeoCompositeShape("sV0AFraR4Sec4","sV0AFraR4bSec4-sV0AChaSec4");
    TGeoCompositeShape *sV0AFraR5Sec4 = new TGeoCompositeShape("sV0AFraR5Sec4","sV0AFraR5bSec4-sV0AChaSec4");
    TGeoVolume *v0AFraR1Sec4 = new TGeoVolume("V0AFraR1Sec4",sV0AFraR1Sec4,medV0AFra);
    TGeoVolume *v0AFraR2Sec4 = new TGeoVolume("V0AFraR2Sec4",sV0AFraR2Sec4,medV0AFra);
    TGeoVolume *v0AFraR3Sec4 = new TGeoVolume("V0AFraR3Sec4",sV0AFraR3Sec4,medV0AFra);
    TGeoVolume *v0AFraR4Sec4 = new TGeoVolume("V0AFraR4Sec4",sV0AFraR4Sec4,medV0AFra);
    TGeoVolume *v0AFraR5Sec4 = new TGeoVolume("V0AFraR5Sec4",sV0AFraR5Sec4,medV0AFra);
    v0AFraR1Sec4->SetLineColor(kV0AColorFra); v0AFraR2Sec4->SetLineColor(kV0AColorFra);
    v0AFraR3Sec4->SetLineColor(kV0AColorFra); v0AFraR4Sec4->SetLineColor(kV0AColorFra);
    v0AFraR5Sec4->SetLineColor(kV0AColorFra);
    v0AFraSec4->AddNode(v0AFraR1Sec4,1); 
    v0AFraSec4->AddNode(v0AFraR2Sec4,1);
    v0AFraSec4->AddNode(v0AFraR3Sec4,1); 
    v0AFraSec4->AddNode(v0AFraR4Sec4,1);
    v0AFraSec4->AddNode(v0AFraR5Sec4,1);
    v0ASec4->AddNode(v0AFraSec4,1);
    
    /// Sensitive scintilator
    TGeoVolume *v0ASciSec4 = new TGeoVolumeAssembly("V0ASciSec4");
    new TGeoTubeSeg( "sV0AR1bSec4", fV0AR0+fV0AFraWd/2.,
		     fV0AR1-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR2bSec4", fV0AR1+fV0AFraWd/2.,
		     fV0AR2-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR3bSec4", fV0AR2+fV0AFraWd/2.,
		     fV0AR3-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    new TGeoTubeSeg( "sV0AR4bSec4", fV0AR3+fV0AFraWd/2.,
		     fV0AR4-fV0AFraWd/2., fV0ASciWd/2., 0, 45);
    TGeoCompositeShape *sV0AR1Sec4 = new TGeoCompositeShape("sV0AR1Sec4","sV0AR1bSec4-sV0AChaSec4");
    TGeoCompositeShape *sV0AR2Sec4 = new TGeoCompositeShape("sV0AR2Sec4","sV0AR2bSec4-sV0AChaSec4");
    TGeoCompositeShape *sV0AR3Sec4 = new TGeoCompositeShape("sV0AR3Sec4","sV0AR3bSec4-sV0AChaSec4");
    TGeoCompositeShape *sV0AR4Sec4 = new TGeoCompositeShape("sV0AR4Sec4","sV0AR4bSec4-sV0AChaSec4");
    TGeoVolume *v0L1Sec4 = new TGeoVolume("V0L1Sec4",sV0AR1Sec4,medV0ASci);
    TGeoVolume *v0L2Sec4 = new TGeoVolume("V0L2Sec4",sV0AR2Sec4,medV0ASci);
    TGeoVolume *v0L3Sec4 = new TGeoVolume("V0L3Sec4",sV0AR3Sec4,medV0ASci);
    TGeoVolume *v0L4Sec4 = new TGeoVolume("V0L4Sec4",sV0AR4Sec4,medV0ASci);
    v0L1Sec4->SetLineColor(kV0AColorSci); v0L2Sec4->SetLineColor(kV0AColorSci);
    v0L3Sec4->SetLineColor(kV0AColorSci); v0L4Sec4->SetLineColor(kV0AColorSci);
    v0ASec4->AddNode(v0L1Sec4,1);
    v0ASec4->AddNode(v0L2Sec4,1);
    v0ASec4->AddNode(v0L3Sec4,1);
    v0ASec4->AddNode(v0L4Sec4,1);      
    
    /// Segment of octagon 
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] =  fV0AR6-fV0AOctH2;           v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = (fV0AR6-fV0AOctH2)*sin45;  v0APts[3+8*i] = (fV0AR6-fV0AOctH2)*sin45;
      v0APts[4+8*i] = fV0AR6*sin45;		 v0APts[5+8*i] = fV0AR6*sin45;
      v0APts[6+8*i] = fV0AR6;			 v0APts[7+8*i] = 0.;
    }
    new TGeoArb8("sV0AOct2Sec4", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoCompositeShape *sV0AOct2FEEBSec4 = new TGeoCompositeShape("sV0AOct2FEEBSec4","sV0AOct2Sec4-sV0AFicFEEBSec4:posFicFEEBSec4-sV0AFicFEEBSec4:posFicFEEBUpSec4-sV0AFicOct2Sec4:posFicOct2Sec4-sV0AFicOct2Sec4:posFicOct2UpSec4");
    TGeoVolume *v0AOct2Sec4 = new TGeoVolume("V0AOct2Sec4", sV0AOct2FEEBSec4,medV0ASup);
    v0AOct2Sec4->SetLineColor(kV0AColorOct);
    TGeoVolume *v0ASupSec4 = new TGeoVolumeAssembly("V0ASupSec4");
    v0ASupSec4->AddNode(v0AOct2Sec4,1);
    v0ASec4->AddNode(v0ASupSec4,1);

    //Bunch of fibers
    v0APts[ 0] = v0APts[ 2] = -13.0;
    v0APts[ 1] = v0APts[ 7] = (fV0ASciWd+fV0AOctWd)/2.-0.01;
    v0APts[ 3] = v0APts[ 5] = (fV0ASciWd+fV0AOctWd)/2.+0.01;
    v0APts[ 4] = v0APts[ 6] = +13.0;
    v0APts[ 8] = v0APts[10] = -10.0;
    v0APts[ 9] = v0APts[15] = 0.;
    v0APts[11] = v0APts[13] = 0.25;
    v0APts[12] = v0APts[14] = +10.0;
    new TGeoArb8("sV0AFib1Sec4", 11.5, v0APts);
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateZ(-90.+22.5);
    TGeoCombiTrans *fibpos1Sec4 = new TGeoCombiTrans("fibpos1Sec4", (fV0AR6-fV0AOctH2+fV0AR5)*cos225/2. - 3.8, (fV0AR6-fV0AOctH2+fV0AR5)*sin225/2.-1.8, 0, rot);
    fibpos1Sec4->RegisterYourself();
    TGeoCompositeShape *sV0AFib1HoleSec4 = new TGeoCompositeShape("sV0AFib1HoleSec4","sV0AFib1Sec4:fibpos1Sec4-sV0AFicR5Sec4"); 
    TGeoVolume *v0AFib1HoleSec4 = new TGeoVolume("V0AFib1HoleSec4",sV0AFib1HoleSec4,medV0AFib);
    v0AFib1HoleSec4->SetLineColor(kV0AColorFib);
    new TGeoArb8("sV0AFib2Sec4", 11.5, v0APts);
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateY(180);
    rot->RotateZ(-90.+22.5);
    TGeoCombiTrans *fibpos2Sec4 = new TGeoCombiTrans("fibpos2Sec4", (fV0AR6-fV0AOctH2+fV0AR5)*cos225/2. - 3.8, (fV0AR6-fV0AOctH2+fV0AR5)*sin225/2.-1.8, 0, rot);
    fibpos2Sec4->RegisterYourself();
    TGeoCompositeShape *sV0AFib2HoleSec4 = new TGeoCompositeShape("sV0AFib2HoleSec4","sV0AFib2Sec4:fibpos2Sec4-sV0AFicR5Sec4");
    TGeoVolume *v0AFib2HoleSec4 = new TGeoVolume("V0AFib2HoleSec4",sV0AFib2HoleSec4,medV0AFib);
    v0AFib2HoleSec4->SetLineColor(kV0AColorFib);
    TGeoVolume *v0AFibSec4 = new TGeoVolumeAssembly("V0AFibSec4");
    v0AFibSec4->AddNode(v0AFib1HoleSec4,1);
    v0AFibSec4->AddNode(v0AFib2HoleSec4,1);
    v0ASec4->AddNode(v0AFibSec4,1); 
    
     /// Plates
    new TGeoTube("sV0ANail1PlaInHoleSec4", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoTube("sV0ANail2PlaInHoleSec4", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoCompositeShape("sV0ANailsPlaInHolesSec4","sV0ANail1PlaInHoleSec4:pos1Sec4+sV0ANail2PlaInHoleSec4:pos2Sec4");
    new TGeoTube("sV0ANail1PlaOuHoleSec4", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoTube("sV0ANail2PlaOuHoleSec4", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoCompositeShape("sV0ANailsPlaOuHolesSec4","sV0ANail1PlaOuHoleSec4:pos1Sec4+sV0ANail2PlaOuHoleSec4:pos2Sec4");
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0;			v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = fV0AR0*sin45;		v0APts[3+8*i] = fV0AR0*sin45;
      v0APts[4+8*i] = fV0AR6 * sin45;	v0APts[5+8*i] = fV0AR6*sin45;
      v0APts[6+8*i] = fV0AR6;		v0APts[7+8*i] = 0.;
    }
    new TGeoArb8("sV0APlaInSec4", (fV0APlaWd-2*fV0APlaAl)/2., v0APts);
    TGeoCompositeShape *sV0APlaInNailsHolesSec4 = new TGeoCompositeShape("sV0APlaInNailsHolesSec4","sV0APlaInSec4-sV0ANailsPlaInHolesSec4-sV0AFicFEEBSec4:posFicFEEBSec4-sV0AFicFEEBSec4:posFicFEEBUpSec4");
    TGeoVolume *v0APlaInNailsHolesSec4 = new TGeoVolume("V0APlaInNailsHolesSec4", sV0APlaInNailsHolesSec4, medV0APlaIn);
    new TGeoArb8("sV0APlaOuSec4", fV0APlaAl/2., v0APts);
    TGeoCompositeShape *sV0APlaOuNailsHolesSec4 = new TGeoCompositeShape("sV0APlaOuNailsHolesSec4","sV0APlaOuSec4-sV0ANailsPlaOuHolesSec4-sV0AFicFEEBSec4:posFicFEEBSec4-sV0AFicFEEBSec4:posFicFEEBUpSec4");  
    TGeoVolume *v0APlaOuNailsHolesSec4 = new TGeoVolume("V0APlaOuNailsHolesSec4", sV0APlaOuNailsHolesSec4, medV0APlaOu);
    v0APlaInNailsHolesSec4->SetLineColor(kV0AColorPlaIn); v0APlaOuNailsHolesSec4->SetLineColor(kV0AColorPlaOu);
    TGeoVolume *v0APlaSec4 = new TGeoVolumeAssembly("V0APlaSec4");
    v0APlaSec4->AddNode(v0APlaInNailsHolesSec4,1);
    v0APlaSec4->AddNode(v0APlaOuNailsHolesSec4,1,new TGeoTranslation(0,0,(fV0APlaWd-fV0APlaAl)/2.));
    v0APlaSec4->AddNode(v0APlaOuNailsHolesSec4,2,new TGeoTranslation(0,0,-(fV0APlaWd-fV0APlaAl)/2.));
    v0ASec4->AddNode(v0APlaSec4,1,new TGeoTranslation(0,0,(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    v0ASec4->AddNode(v0APlaSec4,2,new TGeoTranslation(0,0,-(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    
     /// Non-sensitive scintilator
    new TGeoTubeSeg("sV0AR5S2Sec4", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShapeSec4, 0, 45);
    TGeoCompositeShape *sV0AR5Sec4 = new TGeoCompositeShape("V0AR5Sec4","sV0AR5S2Sec4 - sV0AChaSec4");
    TGeoVolume *v0AR5Sec4 = new TGeoVolume("V0AR5Sec4",sV0AR5Sec4,medV0ASci);
    v0AR5Sec4->SetLineColor(kV0AColorSci);
    v0ASciSec4->AddNode(v0AR5Sec4,1);
    v0ASec4->AddNode(v0ASciSec4,1); 

    /// PMBox
    TGeoVolume* v0APMSec4 = new TGeoVolumeAssembly("V0APMSec4");
    new TGeoBBox("sV0APMB1Sec4", fV0APMBWd/2., fV0APMBHt/2., fV0APMBTh/2.);
    new TGeoBBox("sV0APMB2Sec4", fV0APMBWd/2.-fV0APMBWdW, fV0APMBHt/2.-fV0APMBHtW, fV0APMBTh/2.-fV0APMBThW);
    TGeoCompositeShape *sV0APMBSec4 = new TGeoCompositeShape("sV0APMBSec4","sV0APMB1Sec4-sV0APMB2Sec4");
    TGeoVolume *v0APMBSec4 = new TGeoVolume("V0APMBSec4",sV0APMBSec4, medV0APMAlum);
    v0APMBSec4->SetLineColor(kV0AColorPMA);
    v0APMSec4->AddNode(v0APMBSec4,1);

    /// PMTubes
    TGeoTube *sV0APMT1Sec4 = new TGeoTube("sV0APMT1Sec4", fV0APMTR1, fV0APMTR2, fV0APMTH/2.);
    TGeoVolume *v0APMT1Sec4 = new TGeoVolume("V0APMT1Sec4", sV0APMT1Sec4, medV0APMGlass);
    TGeoTube *sV0APMT2Sec4 = new TGeoTube("sV0APMT2Sec4", fV0APMTR3, fV0APMTR4, fV0APMTH/2.);
    TGeoVolume *v0APMT2Sec4 = new TGeoVolume("V0APMT2Sec4", sV0APMT2Sec4, medV0APMAlum);
    TGeoVolume *v0APMTSec4 = new TGeoVolumeAssembly("V0APMTSec4");
    TGeoTube *sV0APMTTSec4 = new TGeoTube("sV0APMTTSec4", 0., fV0APMTR4, fV0APMTB/2.);
    TGeoVolume *v0APMTTSec4 = new TGeoVolume("V0APMTTSec4", sV0APMTTSec4, medV0APMAlum);
    v0APMT1Sec4->SetLineColor(kV0AColorPMG);
    v0APMT2Sec4->SetLineColor(kV0AColorPMA);
    v0APMTTSec4->SetLineColor(kV0AColorPMA);
    rot = new TGeoRotation("rot", 90, 0, 180, 0, 90, 90);
    v0APMTSec4->AddNode(v0APMT1Sec4,1,rot);
    v0APMTSec4->AddNode(v0APMT2Sec4,1,rot);
    v0APMTSec4->AddNode(v0APMTTSec4,1,new TGeoCombiTrans(0,-(fV0APMTH+fV0APMTB)/2.,0,rot));
    double autoShiftSec4 = (fV0APMBWd-2*fV0APMBWdW)/4.;
    v0APMSec4->AddNode(v0APMTSec4, 1, new TGeoTranslation(-1.5*autoShiftSec4, 0, 0));
    v0APMSec4->AddNode(v0APMTSec4, 2, new TGeoTranslation(-0.5*autoShiftSec4, 0, 0));
    v0APMSec4->AddNode(v0APMTSec4, 3, new TGeoTranslation(+0.5*autoShiftSec4, 0, 0));
    v0APMSec4->AddNode(v0APMTSec4, 4, new TGeoTranslation(+1.5*autoShiftSec4, 0, 0));

    // PM
    rot = new TGeoRotation("rot");
    rot->RotateX(90-fV0APMBAng);
    rot->RotateZ(-90.+22.5);
    double cosAngPMBSec4 = TMath::Cos(fV0APMBAng*TMath::DegToRad());
    double sinAngPMBSec4 = TMath::Sin(fV0APMBAng*TMath::DegToRad());
    double shiftZSec4 = fV0APMBHt/2. * cosAngPMBSec4
      -   ( fV0ASciWd + 2 * fV0AOctWd + 2 * fV0APlaWd )/2.   -   fV0APMBTh/2. * sinAngPMBSec4;
    double shiftRSec4 = fV0AR6  +  fV0AOctH2 + fV0APlaAl;
    v0ASec4->AddNode(v0APMSec4,1, new TGeoCombiTrans( shiftRSec4*cos225+1.07, shiftRSec4*sin225, shiftZSec4, rot));
    
    // Aluminium nails 
    TGeoTube *sV0ANail1Sec4 = new TGeoTube("sV0ANail1Sec4", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail1Sec4 = new TGeoVolume("V0ANail1Sec4", sV0ANail1Sec4, medV0APMAlum);
    v0ANail1Sec4->SetLineColor(kV0AColorPMA);// this is the color for aluminium
    v0ASec4->AddNode(v0ANail1Sec4,1,new TGeoTranslation(42.9, 0.51, 0.0));
    TGeoTube *sV0ANail2Sec4 = new TGeoTube("sV0ANail2Sec4", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail2Sec4 = new TGeoVolume("V0ANail2Sec4", sV0ANail2Sec4, medV0APMAlum);
    v0ANail2Sec4->SetLineColor(kV0AColorPMA);
    v0ASec4->AddNode(v0ANail2Sec4,1,new TGeoTranslation(30.73,29.98,0.0));     
        
    /// Adding sector to v0LE volume
    for(int i=3; i<4; i++) {
       TGeoRotation *rotation = new TGeoRotation("rotation", 90., i*45., 90., 90.+i*45., 0., 0.);
       v0LE->AddNode(v0ASec4,i+1,rotation);  
    }

    //FEEBox
    TGeoVolume* v0AFEE4 = new TGeoVolumeAssembly("V0AFEE4"); 
    v0AFEE4->AddNode(v0AFEEB,1);

    //Mother and daughter boards
    v0AFEE4->AddNode(v0AFEEDaughter, 1, new TGeoTranslation(9.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE4->AddNode(v0AFEEDaughter, 2, new TGeoTranslation(6.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE4->AddNode(v0AFEEDaughter, 3, new TGeoTranslation(3.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE4->AddNode(v0AFEEDaughter, 4, new TGeoTranslation(0.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE4->AddNode(v0AFEEDaughter, 5, new TGeoTranslation(-3.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE4->AddNode(v0AFEEDaughter, 6, new TGeoTranslation(-6.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE4->AddNode(v0AFEEDaughter, 7, new TGeoTranslation(-9.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));
    v0AFEE4->AddNode(v0AFEEDaughter, 8, new TGeoTranslation(-12.0*spacing, fV0AFEEBTh/2.+2.0*fV0APMTB, fV0APMBHtW/2.));   
    v0AFEE4->AddNode(v0AFEEMother, 1, new TGeoTranslation(0.0, 0.0, -fV0AFEEBTh/2.+fV0APMBThW+fV0APMBHtW));
    v0AFEE4->AddNode(v0AFEEHalfMother, 1, new TGeoTranslation(0.0, -(fV0AFEEBTh+fV0APMTB)/2., -2.0*spacing));

    //FEE
    rot = new TGeoRotation("rot");
    rot->RotateX(90);
    rot->RotateZ(-90.0);
    v0LE->AddNode(v0AFEE4,1, new TGeoCombiTrans( -aFEEshiftR2Sec1*cos225-2.0, 0, 7.5, rot));   
    v0LE->AddNode(v0AFEEOct2,4, new TGeoTranslation(-1.0*((aFEEshiftR2Sec1*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0),0.0,0.0));


    //Definition of sector 5
    TGeoVolume *v0ASec5 = new TGeoVolumeAssembly("V0ASec5"); 

    /// For boolean sustraction
    double preShape5 = 0.2;
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = -fV0AR0+fV0AFraWd/2.;  v0APts[1+8*i] = fV0AFraWd/2.;
      v0APts[2+8*i] = -fV0AR0+fV0AFraWd/2.;  v0APts[3+8*i] = -2*preShape5;
      v0APts[4+8*i] = -fV0AR4-fV0AFraWd/2.-preShape5;  v0APts[5+8*i] = -2*preShape5;
      v0APts[6+8*i] = -fV0AR4-fV0AFraWd/2.-preShape5;  v0APts[7+8*i] = fV0AFraWd/2.;
    }
    new TGeoArb8("sV0ACha15",fV0ASciWd/1.5,v0APts);
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = -fV0AR0*cos45+preShape5;
      v0APts[1+8*i] = -(fV0AR0-fV0AFraWd)*sin45+preShape5;
      v0APts[2+8*i] = -(fV0AR0-fV0AFraWd/2.)*cos45+preShape5;
      v0APts[3+8*i] = -(fV0AR0-fV0AFraWd/2.)*sin45;
      v0APts[4+8*i] = -(fV0AR4+fV0AFraWd/2.)*cos45-preShape5;
      v0APts[5+8*i] = -(fV0AR4+fV0AFraWd/2.)*sin45-2.*preShape5;
      v0APts[6+8*i] = -(fV0AR4+fV0AFraWd)*cos45-preShape5;
      v0APts[7+8*i] = -fV0AR4*sin45-preShape5;
    }
    new TGeoArb8("sV0ACha25", fV0ASciWd/2.+2.*preShape5, v0APts);
    new TGeoCompositeShape("sV0ACha125","sV0ACha15+sV0ACha25");
    new TGeoTube("sV0ANail15Hole", 0.0, 0.4, 1.65);
    TGeoTranslation *pos15 = new TGeoTranslation("pos15", -42.9, -0.51, 0.0);
    pos15->RegisterYourself();
    new TGeoTube("sV0ANail25Hole", 0.0, 0.4, 1.65);
    TGeoTranslation *pos25 = new TGeoTranslation("pos25",-30.8,-30.04,0.0); 
    pos25->RegisterYourself();
    TGeoTranslation *pos35 = new TGeoTranslation("pos35",-30.05,-30.79,0.0);  
    pos35->RegisterYourself();
    new TGeoCompositeShape("sV0ANailsHoles5","sV0ANail15Hole:pos15+sV0ANail25Hole:pos25");
    new TGeoCompositeShape("sV0ACha5","sV0ACha125+sV0ANailsHoles5");
    new TGeoTubeSeg("sV0AFicR55", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2. +2*preShape5, 180.0, 225.0);
    new TGeoTube("sV0ANail1PlaInHole5", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoTube("sV0ANail2PlaInHole5", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoTube("sV0ANail3PlaInHole5", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoCompositeShape("sV0ANailsPlaInHoles5","sV0ANail1PlaInHole5:pos15+sV0ANail2PlaInHole5:pos25+sV0ANail3PlaInHole5:pos35");
    new TGeoTube("sV0ANail1PlaOuHole5", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoTube("sV0ANail2PlaOuHole5", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoTube("sV0ANail3PlaOuHole5", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoCompositeShape("sV0ANailsPlaOuHoles5","sV0ANail1PlaOuHole5:pos15+sV0ANail2PlaOuHole5:pos25+sV0ANail3PlaOuHole5:pos35");
    rot = new TGeoRotation("rot");
    rot->RotateX(90);
    rot->RotateZ(-90.0);
    TGeoCombiTrans *posFicFEEBSec5 = new TGeoCombiTrans("posFicFEEBSec5", -aFEEshiftR2Sec1*cos225 - 2.0, 0, 7.5, rot);
    posFicFEEBSec5->RegisterYourself();
    TGeoTranslation *posFicOct2Sec5 = new TGeoTranslation("posFicOct2Sec5",-1.0*((aFEEshiftR2Sec1*cos225 + 2.0) - fV0AFEEBTh/2. - 1.0),0.0,0.0);
    posFicOct2Sec5->RegisterYourself(); 

    /// Frame
    TGeoVolume *v0AFra5 = new TGeoVolumeAssembly("V0AFra5");
    for (int i=0;i<2;i++) { 
      v0APts[0+8*i] = -fV0AR0+fV0AFraWd/2.;  v0APts[1+8*i] = 0.0;
      v0APts[2+8*i] = -fV0AR0+fV0AFraWd/2.;  v0APts[3+8*i] = -fV0AFraWd/8.;
      v0APts[4+8*i] = -fV0AR4-fV0AFraWd/2.;  v0APts[5+8*i] = -fV0AFraWd/8.;
      v0APts[6+8*i] = -fV0AR4-fV0AFraWd/2.;  v0APts[7+8*i] = 0.0;
    }    
    TGeoArb8 *sV0AFraB15 = new TGeoArb8("sV0AFraB15",fV0ASciWd/2.,v0APts);
    TGeoVolume *v0AFraB15 = new TGeoVolume("V0AFraB15",sV0AFraB15,medV0AFra);
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = -fV0AR0*cos45;
      v0APts[1+8*i] = -(fV0AR0-fV0AFraWd)*sin45;
      v0APts[2+8*i] = -(fV0AR0-fV0AFraWd/2.)*cos45;
      v0APts[3+8*i] = -(fV0AR0-fV0AFraWd/2.)*sin45;
      v0APts[4+8*i] = -(fV0AR4+fV0AFraWd/2.)*cos45;
      v0APts[5+8*i] = -(fV0AR4+fV0AFraWd/2.)*sin45;
      v0APts[6+8*i] = -(fV0AR4+fV0AFraWd)*cos45;
      v0APts[7+8*i] = -fV0AR4*sin45;
    }
    TGeoArb8 *sV0AFraB25 = new TGeoArb8("sV0AFraB25", fV0ASciWd/2., v0APts);
    TGeoVolume *v0AFraB25 = new TGeoVolume("V0AFraB25",sV0AFraB25,medV0AFra);
    v0AFraB15->SetLineColor(kV0AColorFra); v0AFraB25->SetLineColor(kV0AColorFra);
    v0AFra5->AddNode(v0AFraB15,1);
    v0AFra5->AddNode(v0AFraB25,1);  // Prefer 2 GeoObjects insted of 3 GeoMovements
    new TGeoTubeSeg( "sV0AFraR1b5", fV0AR0-fV0AFraWd/2.,
		     fV0AR0+fV0AFraWd/2., fV0ASciWd/2., 180.0, 225.0);
    new TGeoTubeSeg( "sV0AFraR2b5", fV0AR1-fV0AFraWd/2.,
		     fV0AR1+fV0AFraWd/2., fV0ASciWd/2., 180.0, 225.0);
    new TGeoTubeSeg( "sV0AFraR3b5", fV0AR2-fV0AFraWd/2.,
		     fV0AR2+fV0AFraWd/2., fV0ASciWd/2., 180.0, 225.0);
    new TGeoTubeSeg( "sV0AFraR4b5", fV0AR3-fV0AFraWd/2.,
		     fV0AR3+fV0AFraWd/2., fV0ASciWd/2., 180.0, 225.0);
    new TGeoTubeSeg( "sV0AFraR5b5", fV0AR4-fV0AFraWd/2.,
		     fV0AR4+fV0AFraWd/2., fV0ASciWd/2., 180.0, 225.0);
    TGeoCompositeShape *sV0AFraR15 = new TGeoCompositeShape("sV0AFraR15","sV0AFraR1b5-sV0ACha5");
    TGeoCompositeShape *sV0AFraR25 = new TGeoCompositeShape("sV0AFraR25","sV0AFraR2b5-sV0ACha5");
    TGeoCompositeShape *sV0AFraR35 = new TGeoCompositeShape("sV0AFraR35","sV0AFraR3b5-sV0ACha5");
    TGeoCompositeShape *sV0AFraR45 = new TGeoCompositeShape("sV0AFraR45","sV0AFraR4b5-sV0ACha5");
    TGeoCompositeShape *sV0AFraR55 = new TGeoCompositeShape("sV0AFraR55","sV0AFraR5b5-sV0ACha5");
    TGeoVolume *v0AFraR15 = new TGeoVolume("V0AFraR15",sV0AFraR15,medV0AFra);
    TGeoVolume *v0AFraR25 = new TGeoVolume("V0AFraR25",sV0AFraR25,medV0AFra);
    TGeoVolume *v0AFraR35 = new TGeoVolume("V0AFraR35",sV0AFraR35,medV0AFra);
    TGeoVolume *v0AFraR45 = new TGeoVolume("V0AFraR45",sV0AFraR45,medV0AFra);
    TGeoVolume *v0AFraR55 = new TGeoVolume("V0AFraR55",sV0AFraR55,medV0AFra);
    v0AFraR15->SetLineColor(kV0AColorFra); v0AFraR25->SetLineColor(kV0AColorFra);
    v0AFraR35->SetLineColor(kV0AColorFra); v0AFraR45->SetLineColor(kV0AColorFra);
    v0AFraR55->SetLineColor(kV0AColorFra);
    v0AFra5->AddNode(v0AFraR15,1);
    v0AFra5->AddNode(v0AFraR25,1);
    v0AFra5->AddNode(v0AFraR35,1);
    v0AFra5->AddNode(v0AFraR45,1);
    v0AFra5->AddNode(v0AFraR55,1);
    v0ASec5->AddNode(v0AFra5,1);

    /// Sensitive scintilator
    TGeoVolume *v0ASci5 = new TGeoVolumeAssembly("V0ASci5");
    new TGeoTubeSeg( "sV0AR1b5", fV0AR0+fV0AFraWd/2.,
		     fV0AR1-fV0AFraWd/2., fV0ASciWd/2., 180.0, 225.0);
    new TGeoTubeSeg( "sV0AR2b5", fV0AR1+fV0AFraWd/2.,
		     fV0AR2-fV0AFraWd/2., fV0ASciWd/2., 180.0, 225.0);
    new TGeoTubeSeg( "sV0AR3b5", fV0AR2+fV0AFraWd/2.,
		     fV0AR3-fV0AFraWd/2., fV0ASciWd/2., 180.0, 225.0);
    new TGeoTubeSeg( "sV0AR4b5", fV0AR3+fV0AFraWd/2.,
		     fV0AR4-fV0AFraWd/2., fV0ASciWd/2., 180.0, 225.0);
    TGeoCompositeShape *sV0AR15 = new TGeoCompositeShape("sV0AR15","sV0AR1b5-sV0ACha5");
    TGeoCompositeShape *sV0AR25 = new TGeoCompositeShape("sV0AR25","sV0AR2b5-sV0ACha5");
    TGeoCompositeShape *sV0AR35 = new TGeoCompositeShape("sV0AR35","sV0AR3b5-sV0ACha5");
    TGeoCompositeShape *sV0AR45 = new TGeoCompositeShape("sV0AR45","sV0AR4b5-sV0ACha5");
    TGeoVolume *v0L15 = new TGeoVolume("V0L15",sV0AR15,medV0ASci);
    TGeoVolume *v0L25 = new TGeoVolume("V0L25",sV0AR25,medV0ASci);
    TGeoVolume *v0L35 = new TGeoVolume("V0L35",sV0AR35,medV0ASci);
    TGeoVolume *v0L45 = new TGeoVolume("V0L45",sV0AR45,medV0ASci);
    v0L15->SetLineColor(kV0AColorSci); v0L25->SetLineColor(kV0AColorSci);
    v0L35->SetLineColor(kV0AColorSci); v0L45->SetLineColor(kV0AColorSci);
    v0ASci5->AddNode(v0L15,1);
    v0ASci5->AddNode(v0L25,1);
    v0ASci5->AddNode(v0L35,1);
    v0ASci5->AddNode(v0L45,1);

     /// Segment of octagon  
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = -fV0AR6+fV0AOctH2;	         v0APts[1+8*i] = 0.;
    v0APts[2+8*i] = -(fV0AR7-fV0AOctH2)*cos654;   v0APts[3+8*i] = -(fV0AR7-fV0AOctH2)*sin654;
    v0APts[4+8*i] = -fV0AR7*cos654;	         v0APts[5+8*i] = -fV0AR7*sin654;
    v0APts[6+8*i] = -fV0AR6;			 v0APts[7+8*i] = 0.;
    }
    new TGeoArb8("sV0AOct25", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoCompositeShape *sV0AOct2FEEB5 = new TGeoCompositeShape("sV0AOct2FEEB5","sV0AOct25-sV0AFicFEEBSec1:posFicFEEBSec5-sV0AFicOct2Sec1:posFicOct2Sec5");
    TGeoVolume *v0AOct25 = new TGeoVolume("V0AOct25", sV0AOct2FEEB5,medV0ASup);
    v0AOct25->SetLineColor(kV0AColorOct);
    TGeoVolume *v0ASup5 = new TGeoVolumeAssembly("V0ASup5");
    v0ASup5->AddNode(v0AOct25,1);
    v0ASec5->AddNode(v0ASup5,1);

    //Bunch of fibers
    v0APts[ 0] = v0APts[ 2] = -14.0;
    v0APts[ 1] = v0APts[ 7] = (fV0ASciWd+fV0AOctWd)/2.-0.01;
    v0APts[ 3] = v0APts[ 5] = (fV0ASciWd+fV0AOctWd)/2.+0.01;
    v0APts[ 4] = v0APts[ 6] = +14.0;
    v0APts[ 8] = v0APts[10] = -10.0;
    v0APts[ 9] = v0APts[15] = 0.;
    v0APts[11] = v0APts[13] = 0.25;
    v0APts[12] = v0APts[14] = +10.0;
    new TGeoArb8("sV0AFib15", 11.8, v0APts); 
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateZ(90+22.5);
    TGeoCombiTrans *fib15pos = new TGeoCombiTrans("fib15pos", -(fV0AR6-fV0AOctH2+fV0AR5)*cos225/2. +
    3.3, -(fV0AR6-fV0AOctH2+fV0AR5)*sin225/2. + 1.5, 0, rot);
    fib15pos->RegisterYourself();
    TGeoCompositeShape *sV0AFib15Hole = new TGeoCompositeShape("sV0AFib15Hole", "sV0AFib15:fib15pos-sV0AFicR55");
    TGeoVolume *v0AFib15Hole = new TGeoVolume("V0AFib15",sV0AFib15Hole,medV0AFib);
    v0AFib15Hole->SetLineColor(kV0AColorFib);
    new TGeoArb8("sV0AFib25", 11.8, v0APts);
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateY(180);
    rot->RotateZ(90+22.5);
    TGeoCombiTrans *fib25pos = new TGeoCombiTrans("fib25pos", -(fV0AR6-fV0AOctH2+fV0AR5)*cos225/2. +
    3.3, -(fV0AR6-fV0AOctH2+fV0AR5)*sin225/2. + 1.5, 0, rot);
    fib25pos->RegisterYourself();
    TGeoCompositeShape *sV0AFib25Hole = new TGeoCompositeShape("sV0AFib25Hole","sV0AFib25:fib25pos-sV0AFicR55");
    TGeoVolume *v0AFib25Hole = new TGeoVolume("V0AFib25Hole",sV0AFib25Hole,medV0AFib);
    v0AFib25Hole->SetLineColor(kV0AColorFib);
    TGeoVolume *v0AFib5 = new TGeoVolumeAssembly("V0AFib5");    
    v0AFib5->AddNode(v0AFib15Hole,1);
    v0AFib5->AddNode(v0AFib25Hole,1);
    v0ASec5->AddNode(v0AFib5,1);
            
    /// Non-sensitive scintilator
    new TGeoTubeSeg("sV0AR5S25", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShape5, 180.0, 225.0);
    TGeoCompositeShape *sV0AR55 = new TGeoCompositeShape("V0AR55","sV0AR5S25 - sV0ACha5");
    TGeoVolume *v0AR55 = new TGeoVolume("V0AR55",sV0AR55,medV0ASci);
    v0AR55->SetLineColor(kV0AColorSci);
    v0ASci5->AddNode(v0AR55,1);
    v0ASec5->AddNode(v0ASci5,1);

    /// Plates 
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = -fV0AR0;			v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = -fV0AR0*cos654;		v0APts[3+8*i] = -fV0AR0*sin654;
      v0APts[4+8*i] = -fV0AR7*cos654;    	v0APts[5+8*i] = -fV0AR7*sin654;
      v0APts[6+8*i] = -fV0AR6;    		v0APts[7+8*i] = 0.;
    }
    new TGeoArb8("sV0APlaIn5", (fV0APlaWd-2*fV0APlaAl)/2., v0APts);
    TGeoCompositeShape *sV0APlaInNailsHoles5 = new TGeoCompositeShape("sV0APlaInNailsHoles5","sV0APlaIn5-sV0ANailsPlaInHoles5-sV0AFicFEEBSec1:posFicFEEBSec5");
    TGeoVolume *v0APlaInNailsHoles5 = new TGeoVolume("V0APlaInNailsHoles5", sV0APlaInNailsHoles5, medV0APlaIn);
    new TGeoArb8("sV0APlaOu5", fV0APlaAl/2., v0APts);
    TGeoCompositeShape *sV0APlaOuNailsHoles5 = new TGeoCompositeShape("sV0APlaOuNailsHoles5","sV0APlaOu5-sV0ANailsPlaOuHoles5-sV0AFicFEEBSec1:posFicFEEBSec5");  
    TGeoVolume *v0APlaOuNailsHoles5 = new TGeoVolume("V0APlaOuNailsHoles5", sV0APlaOuNailsHoles5, medV0APlaOu);
    v0APlaInNailsHoles5->SetLineColor(kV0AColorPlaIn); v0APlaOuNailsHoles5->SetLineColor(kV0AColorPlaOu);
    TGeoVolume *v0APla5 = new TGeoVolumeAssembly("V0APla5");
    v0APla5->AddNode(v0APlaInNailsHoles5,1);
    v0APla5->AddNode(v0APlaOuNailsHoles5,1,new TGeoTranslation(0,0,(fV0APlaWd-fV0APlaAl)/2.));
    v0APla5->AddNode(v0APlaOuNailsHoles5,2,new TGeoTranslation(0,0,-(fV0APlaWd-fV0APlaAl)/2.));
    v0ASec5->AddNode(v0APla5,1,new TGeoTranslation(0,0,(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    v0ASec5->AddNode(v0APla5,2,new TGeoTranslation(0,0,-(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    
    /// PMBox 
    TGeoVolume* v0APM5 = new TGeoVolumeAssembly("V0APM5");
    new TGeoBBox("sV0APMB15", fV0APMBWd/2., fV0APMBHt/2., fV0APMBTh/2.);
    new TGeoBBox("sV0APMB25", fV0APMBWd/2.-fV0APMBWdW, fV0APMBHt/2.-fV0APMBHtW, fV0APMBTh/2.-fV0APMBThW);
    TGeoCompositeShape *sV0APMB5 = new TGeoCompositeShape("sV0APMB5","sV0APMB15-sV0APMB25");
    TGeoVolume *v0APMB5 = new TGeoVolume("V0APMB5",sV0APMB5, medV0APMAlum);
    v0APMB5->SetLineColor(kV0AColorPMA);
    v0APM5->AddNode(v0APMB5,1);

    /// PMTubes 
    TGeoTube *sV0APMT15 = new TGeoTube("sV0APMT15", fV0APMTR1, fV0APMTR2, fV0APMTH/2.);
    TGeoVolume *v0APMT15 = new TGeoVolume("V0APMT15", sV0APMT15, medV0APMGlass);
    TGeoTube *sV0APMT25 = new TGeoTube("sV0APMT25", fV0APMTR3, fV0APMTR4, fV0APMTH/2.);
    TGeoVolume *v0APMT25 = new TGeoVolume("V0APMT25", sV0APMT25, medV0APMAlum);
    TGeoVolume *v0APMT5 = new TGeoVolumeAssembly("V0APMT5");
    TGeoTube *sV0APMTT5 = new TGeoTube("sV0APMTT5", 0., fV0APMTR4, fV0APMTB/2.);
    TGeoVolume *v0APMTT5 = new TGeoVolume("V0APMTT5", sV0APMTT5, medV0APMAlum);
    v0APMT5->SetLineColor(kV0AColorPMG);
    v0APMT25->SetLineColor(kV0AColorPMA);
    v0APMTT5->SetLineColor(kV0AColorPMA);
    rot = new TGeoRotation("rot", 90, 0, 180, 0, 90, 90);
    v0APMT5->AddNode(v0APMT15,1,rot);
    v0APMT5->AddNode(v0APMT25,1,rot);
    v0APMT5->AddNode(v0APMTT5,1,new TGeoCombiTrans(0,(fV0APMTH+fV0APMTB)/2.,0,rot));
    double autoShift5 = (fV0APMBWd-2*fV0APMBWdW)/4.;
    v0APM5->AddNode(v0APMT5, 1, new TGeoTranslation(-1.5*autoShift5, 0, 0));
    v0APM5->AddNode(v0APMT5, 2, new TGeoTranslation(-0.5*autoShift5, 0, 0));
    v0APM5->AddNode(v0APMT5, 3, new TGeoTranslation(+0.5*autoShift5, 0, 0));
    v0APM5->AddNode(v0APMT5, 4, new TGeoTranslation(+1.5*autoShift5, 0, 0));

    /// PM 
    rot = new TGeoRotation("rot");
    rot->RotateX(-90+30);
    rot->RotateY(0); 
    rot->RotateZ(-65-3);
    double cosAngPMB5 = TMath::Cos(fV0APMBAng*TMath::DegToRad());
    double sinAngPMB5 = TMath::Sin(fV0APMBAng*TMath::DegToRad());
    double shiftZ5 = fV0APMBHt/2. * cosAngPMB5
      -   ( fV0ASciWd + 2 * fV0AOctWd + 2 * fV0APlaWd )/2.   -   fV0APMBTh/2. * sinAngPMB5;
    double shiftR5 = fV0AR6  +  fV0AOctH2 + fV0APlaAl;
    v0ASec5->AddNode(v0APM5,1, new TGeoCombiTrans( -shiftR5*cos225-1.3, -shiftR5*sin225, shiftZ5, rot));

    // Aluminium nails
    TGeoTube *sV0ANail51 = new TGeoTube("sV0ANail51", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail51 = new TGeoVolume("V0ANail51", sV0ANail51, medV0APMAlum);
    v0ANail51->SetLineColor(kV0AColorPMA);// this is the color for aluminium
    v0ASec5->AddNode(v0ANail51,1,new TGeoTranslation(-42.9,-0.51,0.0));
    TGeoTube *sV0ANail52 = new TGeoTube("sV0ANail52", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail52 = new TGeoVolume("V0ANail52", sV0ANail52, medV0APMAlum);
    v0ANail52->SetLineColor(kV0AColorPMA);
    v0ASec5->AddNode(v0ANail52,1,new TGeoTranslation(-30.8,-30.04,0.0)); 
            
    // Adding sector to v0LE volume
    v0LE->AddNode(v0ASec5, 1);  
    

    //Definition of  sector 6
    TGeoVolume *v0ASec6 = new TGeoVolumeAssembly("V0ASec6");

    /// For boolean sustraction
    double preShape6 = 0.2;
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = -preShape6;                      v0APts[1+8*i] = -fV0AR0+fV0AFraWd/2.-preShape6;
    v0APts[2+8*i] = 0.0;                    v0APts[3+8*i] = -fV0AR0+fV0AFraWd/2.-preShape6;
    v0APts[4+8*i] = 0.0;                    v0APts[5+8*i] = -fV0AR4-fV0AFraWd/2.+preShape6;
    v0APts[6+8*i] = -preShape6;                      v0APts[7+8*i] = -fV0AR4-fV0AFraWd/2.+preShape6;
    }
    new TGeoArb8("sV0ACha16",fV0ASciWd/1.5,v0APts);
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = -fV0AR0*cos45+preShape6;
    v0APts[1+8*i] = -(fV0AR0-fV0AFraWd)*sin45+preShape6;
    v0APts[2+8*i] = -(fV0AR0-fV0AFraWd/2.)*cos45+preShape6;
    v0APts[3+8*i] = -(fV0AR0-fV0AFraWd/2.)*sin45;
    v0APts[4+8*i] = -(fV0AR4+fV0AFraWd/2.)*cos45-preShape6;
    v0APts[5+8*i] = -(fV0AR4+fV0AFraWd/2.)*sin45-preShape6;
    v0APts[6+8*i] = -(fV0AR4+fV0AFraWd)*cos45-preShape6;
    v0APts[7+8*i] = -fV0AR4*sin45-preShape6;
    }
    new TGeoArb8("sV0ACha26", fV0ASciWd/2.+2.*preShape6, v0APts);
    new TGeoCompositeShape("sV0ACha126","sV0ACha16+sV0ACha26");
    new TGeoTube("sV0ANail16Hole", 0.0, 0.4, 1.65);
    TGeoTranslation *pos16 = new TGeoTranslation("pos16",-0.51,-42.9,0.0);
    pos16->RegisterYourself();
    new TGeoTube("sV0ANail26Hole", 0.0, 0.4, 1.65);
    TGeoTranslation *pos26 = new TGeoTranslation("pos26",-30.05,-30.79,0.0);  
    pos26->RegisterYourself();
    new TGeoCompositeShape("sV0ANailsHoles6","sV0ANail16Hole:pos16+sV0ANail26Hole:pos26");
    new TGeoCompositeShape("sV0ACha6","sV0ACha126+sV0ANailsHoles6");
    new TGeoTubeSeg("sV0AFicR56", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShape6, 225, 270.0);
    new TGeoTube("sV0ANail1PlaInHole6", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);    
    new TGeoTube("sV0ANail1PlaOuHole6", 0.0, 0.4, (fV0APlaAl)/2.);
      
    /// Frame
    TGeoVolume *v0AFra6 = new TGeoVolumeAssembly("V0AFra6");
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = -fV0AFraWd/2.;         v0APts[1+8*i] = -fV0AR0-fV0AFraWd/2.;
    v0APts[2+8*i] = 0.;                    v0APts[3+8*i] = -fV0AR0-fV0AFraWd/2.;
    v0APts[4+8*i] = 0.;                    v0APts[5+8*i] = -fV0AR4+fV0AFraWd/2.;
    v0APts[6+8*i] = -fV0AFraWd/2.;         v0APts[7+8*i] = -fV0AR4+fV0AFraWd/2.;
    }
    TGeoArb8 *sV0AFraB16 = new TGeoArb8("sV0AFraB16",fV0ASciWd/2.,v0APts);
    TGeoVolume *v0AFraB16 = new TGeoVolume("V0AFraB16",sV0AFraB16,medV0AFra);
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = -(fV0AR0+fV0AFraWd/2.)*cos45;
    v0APts[1+8*i] = -(fV0AR0+fV0AFraWd/2.)*sin45;
    v0APts[2+8*i] = -fV0AR0*cos45;
    v0APts[3+8*i] = -(fV0AR0+fV0AFraWd)*sin45;
    v0APts[4+8*i] = -(fV0AR4-fV0AFraWd)*cos45;
    v0APts[5+8*i] = -fV0AR4*sin45;
    v0APts[6+8*i] = -(fV0AR4-fV0AFraWd/6.)*cos45;
    v0APts[7+8*i] = -(fV0AR4-fV0AFraWd/2.)*sin45;
    }
    TGeoArb8 *sV0AFraB26 = new TGeoArb8("sV0AFraB26", fV0ASciWd/2., v0APts);
    TGeoVolume *v0AFraB26 = new TGeoVolume("V0AFraB26",sV0AFraB26,medV0AFra);
    v0AFraB16->SetLineColor(kV0AColorFra); v0AFraB26->SetLineColor(kV0AColorFra);
    v0AFra6->AddNode(v0AFraB16,1);
    v0AFra6->AddNode(v0AFraB26,1);  // Prefer 2 GeoObjects insted of 3 GeoMovements
    new TGeoTubeSeg( "sV0AFraR1b6", fV0AR0-fV0AFraWd/2.,
    	     fV0AR0+fV0AFraWd/2., fV0ASciWd/2., 225.0, 270.0);
    new TGeoTubeSeg( "sV0AFraR2b6", fV0AR1-fV0AFraWd/2.,
    	     fV0AR1+fV0AFraWd/2., fV0ASciWd/2., 225.0, 270.0);
    new TGeoTubeSeg( "sV0AFraR3b6", fV0AR2-fV0AFraWd/2.,
    	     fV0AR2+fV0AFraWd/2., fV0ASciWd/2., 225.0, 270.0);
    new TGeoTubeSeg( "sV0AFraR4b6", fV0AR3-fV0AFraWd/2.,
    	     fV0AR3+fV0AFraWd/2., fV0ASciWd/2., 225.0, 270.0);
    new TGeoTubeSeg( "sV0AFraR5b6", fV0AR4-fV0AFraWd/2.,
    	     fV0AR4+fV0AFraWd/2., fV0ASciWd/2., 225.0, 270.0);
    TGeoCompositeShape *sV0AFraR16 = new TGeoCompositeShape("sV0AFraR16","sV0AFraR1b6-sV0ACha6");
    TGeoCompositeShape *sV0AFraR26 = new TGeoCompositeShape("sV0AFraR26","sV0AFraR2b6-sV0ACha6");
    TGeoCompositeShape *sV0AFraR36 = new TGeoCompositeShape("sV0AFraR36","sV0AFraR3b6-sV0ACha6");
    TGeoCompositeShape *sV0AFraR46 = new TGeoCompositeShape("sV0AFraR46","sV0AFraR4b6-sV0ACha6");
    TGeoCompositeShape *sV0AFraR56 = new TGeoCompositeShape("sV0AFraR56","sV0AFraR5b6-sV0ACha6");
    TGeoVolume *v0AFraR16 = new TGeoVolume("V0AFraR16",sV0AFraR16,medV0AFra);
    TGeoVolume *v0AFraR26 = new TGeoVolume("V0AFraR26",sV0AFraR26,medV0AFra);
    TGeoVolume *v0AFraR36 = new TGeoVolume("V0AFraR36",sV0AFraR36,medV0AFra);
    TGeoVolume *v0AFraR46 = new TGeoVolume("V0AFraR46",sV0AFraR46,medV0AFra);
    TGeoVolume *v0AFraR56 = new TGeoVolume("V0AFraR56",sV0AFraR56,medV0AFra);
    v0AFraR16->SetLineColor(kV0AColorFra); v0AFraR26->SetLineColor(kV0AColorFra);
    v0AFraR36->SetLineColor(kV0AColorFra); v0AFraR46->SetLineColor(kV0AColorFra);
    v0AFraR56->SetLineColor(kV0AColorFra);
    v0AFra6->AddNode(v0AFraR16,1);
    v0AFra6->AddNode(v0AFraR26,1);
    v0AFra6->AddNode(v0AFraR36,1);
    v0AFra6->AddNode(v0AFraR46,1);
    v0AFra6->AddNode(v0AFraR56,1);
    v0ASec6->AddNode(v0AFra6,1);

    /// Sensitive scintilator
    TGeoVolume *v0ASci6 = new TGeoVolumeAssembly("V0ASci6");
    new TGeoTubeSeg( "sV0AR1b6", fV0AR0+fV0AFraWd/2.,
    		     fV0AR1-fV0AFraWd/2., fV0ASciWd/2., 225.0, 270.0);
    new TGeoTubeSeg( "sV0AR2b6", fV0AR1+fV0AFraWd/2.,
    	     fV0AR2-fV0AFraWd/2., fV0ASciWd/2., 225.0, 270.0);
    new TGeoTubeSeg( "sV0AR3b6", fV0AR2+fV0AFraWd/2.,
    	     fV0AR3-fV0AFraWd/2., fV0ASciWd/2., 225.0, 270.0);
    new TGeoTubeSeg( "sV0AR4b6", fV0AR3+fV0AFraWd/2.,
    	     fV0AR4-fV0AFraWd/2., fV0ASciWd/2., 225.0, 270.0);
    TGeoCompositeShape *sV0AR16 = new TGeoCompositeShape("sV0AR16","sV0AR1b6-sV0ACha6");
    TGeoCompositeShape *sV0AR26 = new TGeoCompositeShape("sV0AR26","sV0AR2b6-sV0ACha6");
    TGeoCompositeShape *sV0AR36 = new TGeoCompositeShape("sV0AR36","sV0AR3b6-sV0ACha6");
    TGeoCompositeShape *sV0AR46 = new TGeoCompositeShape("sV0AR46","sV0AR4b6-sV0ACha6");
    TGeoVolume *v0L16 = new TGeoVolume("V0L16",sV0AR16,medV0ASci);
    TGeoVolume *v0L26 = new TGeoVolume("V0L26",sV0AR26,medV0ASci);
    TGeoVolume *v0L36 = new TGeoVolume("V0L36",sV0AR36,medV0ASci);
    TGeoVolume *v0L46 = new TGeoVolume("V0L46",sV0AR46,medV0ASci);
    v0L16->SetLineColor(kV0AColorSci); v0L26->SetLineColor(kV0AColorSci);
    v0L36->SetLineColor(kV0AColorSci); v0L46->SetLineColor(kV0AColorSci);
    v0ASci6->AddNode(v0L16,1);
    v0ASci6->AddNode(v0L26,1);
    v0ASci6->AddNode(v0L36,1);
    v0ASci6->AddNode(v0L46,1);
    
    // Bunch of fibers
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = -10.0;	    v0APts[1+8*i] = 13.1;  
    v0APts[2+8*i] = 10.0;           v0APts[3+8*i] = 13.1;   
    v0APts[4+8*i] = 8.0;            v0APts[5+8*i] = -29.0;  
    v0APts[6+8*i] = -12.0;	    v0APts[7+8*i] = -12.0;  
    }   
    new TGeoArb8("sV0AFib16", 0.01, v0APts);      
    rot = new TGeoRotation("rot");
    rot->RotateX(2.0); 
    rot->RotateY(180.0);
    rot->RotateZ(90+22.5);
    TGeoCombiTrans *fib16pos = new TGeoCombiTrans("fib16pos", -40.0 + 3.3, -50.0 + 1.5, 0.5, rot);
    fib16pos->RegisterYourself();
    TGeoCompositeShape *sV0AFib16Hole = new TGeoCompositeShape("sV0AFib16Hole", "sV0AFib16:fib16pos-sV0AFicR56");
    TGeoVolume *v0AFib16Hole = new TGeoVolume("V0AFib16Hole",sV0AFib16Hole,medV0AFib);
    v0AFib16Hole->SetLineColor(kV0AColorFib);
    new TGeoArb8("sV0AFib26", 0.01, v0APts);      
    rot = new TGeoRotation("rot");
    rot->RotateX(-2.0); 
    rot->RotateY(180.0); 
    rot->RotateZ(90+22.5);
    TGeoCombiTrans *fib26pos = new TGeoCombiTrans("fib26pos", -40.0 + 3.3, -50.0 + 1.5, -0.5, rot);
    fib26pos->RegisterYourself();
    TGeoCompositeShape *sV0AFib26Hole = new TGeoCompositeShape("sV0AFib26Hole", "sV0AFib26:fib26pos-sV0AFicR56");
    TGeoVolume *v0AFib26Hole = new TGeoVolume("V0AFib26Hole",sV0AFib26Hole,medV0AFib);
    v0AFib26Hole->SetLineColor(kV0AColorFib);
    TGeoVolume *v0AFib6 = new TGeoVolumeAssembly("V0AFib6");
    v0AFib6->AddNode(v0AFib16Hole,1); 
    v0AFib6->AddNode(v0AFib26Hole,1);
    v0ASec6->AddNode(v0AFib6,1);

    /// Non-sensitive scintilator
    new TGeoTubeSeg("sV0AR5S26", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShape6, 225.0, 270.0);
    TGeoCompositeShape *sV0AR56 = new TGeoCompositeShape("V0AR56","sV0AR5S26 - sV0ACha6");
    TGeoVolume *v0AR56 = new TGeoVolume("V0AR56",sV0AR56,medV0ASci);
    v0AR56->SetLineColor(kV0AColorSci);
    v0ASci6->AddNode(v0AR56,1);
    v0ASec6->AddNode(v0ASci6,1);

    /// Segment of octagon   
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = 0.;                 	 v0APts[1+8*i] = -(fV0AR7-fV0AOctH2)*sin654;
    v0APts[2+8*i] = 0.;                          v0APts[3+8*i] = -fV0AR7*sin654;
    v0APts[4+8*i] = -fV0AR7*cos654;     	 v0APts[5+8*i] = -fV0AR7*sin654;
    v0APts[6+8*i] = -(fV0AR7-fV0AOctH2)*cos654;	 v0APts[7+8*i] = -(fV0AR7-fV0AOctH2)*sin654;
    }
    TGeoArb8 *sV0AOct26 = new TGeoArb8("sV0AOct26", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoVolume *v0AOct26 = new TGeoVolume("V0AOct26", sV0AOct26,medV0ASup);
    v0AOct26->SetLineColor(kV0AColorOct);
    TGeoVolume *v0ASup6 = new TGeoVolumeAssembly("V0ASup6");
    v0ASup6->AddNode(v0AOct26,1);
    v0ASec6->AddNode(v0ASup6,1);

    /// Plates
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = 0.;       		v0APts[1+8*i] = -fV0AR0;
    v0APts[2+8*i] = 0.;        		v0APts[3+8*i] = -fV0AR7*sin654;
    v0APts[4+8*i] = -fV0AR7*cos654;    	v0APts[5+8*i] = -fV0AR7*sin654;
    v0APts[6+8*i] = -fV0AR0*cos654;  		v0APts[7+8*i] = -fV0AR0*sin654;
    }
    new TGeoArb8("sV0APlaIn6", (fV0APlaWd-2*fV0APlaAl)/2., v0APts);
    TGeoCompositeShape *sV0APlaInNailsHoles6 = new TGeoCompositeShape("sV0APlaInNailsHoles6","sV0APlaIn6-sV0ANail1PlaInHole6:pos16");
    TGeoVolume *v0APlaInNailsHoles6 = new TGeoVolume("V0APlaInNailsHoles6", sV0APlaInNailsHoles6, medV0APlaIn);
    new TGeoArb8("sV0APlaOu6", fV0APlaAl/2., v0APts);
    TGeoCompositeShape *sV0APlaOuNailsHoles6 = new TGeoCompositeShape("sV0APlaOuNailsHoles6","sV0APlaOu6-sV0ANail1PlaOuHole6:pos16"); 
    TGeoVolume *v0APlaOuNailsHoles6 = new TGeoVolume("V0APlaOuNailsHoles6", sV0APlaOuNailsHoles6, medV0APlaOu);
    v0APlaInNailsHoles6->SetLineColor(kV0AColorPlaIn); v0APlaOuNailsHoles6->SetLineColor(kV0AColorPlaOu);
    TGeoVolume *v0APla6 = new TGeoVolumeAssembly("V0APla6");
    v0APla6->AddNode(v0APlaInNailsHoles6,1);
    v0APla6->AddNode(v0APlaOuNailsHoles6,1,new TGeoTranslation(0,0,(fV0APlaWd-fV0APlaAl)/2.));
    v0APla6->AddNode(v0APlaOuNailsHoles6,2,new TGeoTranslation(0,0,-(fV0APlaWd-fV0APlaAl)/2.));
    v0ASec6->AddNode(v0APla6,1,new TGeoTranslation(0,0,(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    v0ASec6->AddNode(v0APla6,2,new TGeoTranslation(0,0,-(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    
    /// PMBox  
    TGeoVolume* v0APM6 = new TGeoVolumeAssembly("V0APM6");
    new TGeoBBox("sV0APMB16", fV0APMBWd/2., fV0APMBHt/2., fV0APMBTh/2.);
    new TGeoBBox("sV0APMB26", fV0APMBWd/2.-fV0APMBWdW, fV0APMBHt/2.-fV0APMBHtW, fV0APMBTh/2.-fV0APMBThW);
    TGeoCompositeShape *sV0APMB6 = new TGeoCompositeShape("sV0APMB6","sV0APMB16-sV0APMB26");
    TGeoVolume *v0APMB6 = new TGeoVolume("V0APMB6",sV0APMB6, medV0APMAlum);
    v0APMB6->SetLineColor(kV0AColorPMA);
    v0APM6->AddNode(v0APMB6,1);

    /// PMTubes 
    TGeoTube *sV0APMT16 = new TGeoTube("sV0APMT16", fV0APMTR1, fV0APMTR2, fV0APMTH/2.);
    TGeoVolume *v0APMT16 = new TGeoVolume("V0APMT16", sV0APMT16, medV0APMGlass);
    TGeoTube *sV0APMT26 = new TGeoTube("sV0APMT26", fV0APMTR3, fV0APMTR4, fV0APMTH/2.);
    TGeoVolume *v0APMT26 = new TGeoVolume("V0APMT26", sV0APMT26, medV0APMAlum);
    TGeoVolume *v0APMT6 = new TGeoVolumeAssembly("V0APMT6"); 
    TGeoTube *sV0APMTT6 = new TGeoTube("sV0APMTT6", 0., fV0APMTR4, fV0APMTB/2.);
    TGeoVolume *v0APMTT6 = new TGeoVolume("V0APMTT6", sV0APMTT6, medV0APMAlum);
    v0APMT6->SetLineColor(kV0AColorPMG);
    v0APMT26->SetLineColor(kV0AColorPMA);
    v0APMTT6->SetLineColor(kV0AColorPMA);
    rot = new TGeoRotation("rot", 90, 0, 180, 0, 90, 90);
    v0APMT6->AddNode(v0APMT16,1,rot);
    v0APMT6->AddNode(v0APMT26,1,rot);
    v0APMT6->AddNode(v0APMTT6,1,new TGeoCombiTrans(0,(fV0APMTH+fV0APMTB)/2.,0,rot));
    double autoShift6 = (fV0APMBWd-2*fV0APMBWdW)/4.;
    v0APM6->AddNode(v0APMT6, 1, new TGeoTranslation(-1.5*autoShift6, 0, 0));
    v0APM6->AddNode(v0APMT6, 2, new TGeoTranslation(-0.5*autoShift6, 0, 0));
    v0APM6->AddNode(v0APMT6, 3, new TGeoTranslation(+0.5*autoShift6, 0, 0));
    v0APM6->AddNode(v0APMT6, 4, new TGeoTranslation(+1.5*autoShift6, 0, 0));

    /// PM 
    rot = new TGeoRotation("rot");
    rot->RotateX(-90+30);
    rot->RotateY(0);
    rot->RotateZ(-65-3);
    double cosAngPMB6 = TMath::Cos(fV0APMBAng*TMath::DegToRad());
    double sinAngPMB6 = TMath::Sin(fV0APMBAng*TMath::DegToRad());
    double shiftZ6 = fV0APMBHt/2. * cosAngPMB6
      -   ( fV0ASciWd + 2 * fV0AOctWd + 2 * fV0APlaWd )/2.   -   fV0APMBTh/2. * sinAngPMB6;
    double shiftR6 = fV0AR6  + fV0AR1 + fV0AOctWd + fV0APlaAl/3.;
    v0ASec6->AddNode(v0APM6,1, new TGeoCombiTrans( -shiftR6*cos45-1.3, -shiftR6*sin45, shiftZ6, rot));   
    
    /// Support
    TGeoBBox *sV0ASuppbl = new TGeoBBox("sV0ASuppbl", 2.0, 18.13, fV0ASciWd/2.);       
    TGeoVolume *v0ASuppbl = new TGeoVolume("V0ASuppbl", sV0ASuppbl, medV0ASup);
    v0ASuppbl->SetLineColor(kV0AColorOct);
    v0ASec6->AddNode(v0ASuppbl,1,new TGeoTranslation(-2.0,-63.64,0.0));
    
    // Aluminium nail
    TGeoTube *sV0ANail61 = new TGeoTube("sV0ANail61", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail61 = new TGeoVolume("V0ANail61", sV0ANail61, medV0APMAlum);
    v0ANail61->SetLineColor(kV0AColorPMA);// this is the color for aluminium
    v0ASec6->AddNode(v0ANail61,1,new TGeoTranslation(-0.51,-42.9,0.0));  
    TGeoTube *sV0ANail62 = new TGeoTube("sV0ANail62", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail62 = new TGeoVolume("V0ANail62", sV0ANail62, medV0APMAlum);
    v0ANail62->SetLineColor(kV0AColorPMA);
    v0ASec6->AddNode(v0ANail62,1,new TGeoTranslation(-30.05,-30.79,0.0));  

    // Adding sector to v0LE volume
    v0LE->AddNode(v0ASec6, 1); 
 
     
     //Definition of sector 7
    TGeoVolume *v0ASec7 = new TGeoVolumeAssembly("V0ASec7");

    /// For boolean sustraction
    double preShape7 = 0.2;
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = 0.0;                      v0APts[1+8*i] = -fV0AR0+fV0AFraWd/2.-preShape7;
    v0APts[2+8*i] = fV0AFraWd/2.;                    v0APts[3+8*i] = -fV0AR0+fV0AFraWd/2.-preShape7;
    v0APts[4+8*i] = fV0AFraWd/2.;                    v0APts[5+8*i] = -fV0AR4-fV0AFraWd/2.+preShape7;
    v0APts[6+8*i] = 0.0;                      v0APts[7+8*i] = -fV0AR4-fV0AFraWd/2.+preShape7;
    }
    new TGeoArb8("sV0ACha17",fV0ASciWd/1.5,v0APts);
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = (fV0AR0-fV0AFraWd/2.)*cos45-preShape7;
    v0APts[1+8*i] =  -(fV0AR0-fV0AFraWd/2.)*sin45;
    v0APts[2+8*i] = fV0AR0*cos45-preShape7;
    v0APts[3+8*i] =  -(fV0AR0-fV0AFraWd)*sin45+preShape7;
    v0APts[4+8*i] = (fV0AR4+fV0AFraWd)*cos45+preShape7;
    v0APts[5+8*i] =  -fV0AR4*sin45-preShape7;
    v0APts[6+8*i] = (fV0AR4+fV0AFraWd/2.)*cos45+preShape7;
    v0APts[7+8*i] =  -(fV0AR4+fV0AFraWd/2.)*sin45-2.*preShape7;
    }
    new TGeoArb8("sV0ACha27", fV0ASciWd/2.+2.*preShape7, v0APts);
    new TGeoCompositeShape("sV0ACha127","sV0ACha17+sV0ACha27");
    new TGeoTube("sV0ANail17Hole", 0.0, 0.4, 1.65);
    TGeoTranslation *pos17 = new TGeoTranslation("pos17",0.51,-42.9,0.0);
    pos17->RegisterYourself();
    new TGeoTube("sV0ANail27Hole", 0.0, 0.4, 1.65);
    TGeoTranslation *pos27 = new TGeoTranslation("pos27",30.05,-30.79,0.0);   
    pos27->RegisterYourself();
    new TGeoCompositeShape("sV0ANailsHoles7","sV0ANail17Hole:pos17+sV0ANail27Hole:pos27");
    new TGeoCompositeShape("sV0ACha7","sV0ACha127+sV0ANailsHoles7");
    new TGeoTubeSeg("sV0AFicR57", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShape7, 270.0, 315.0);
    new TGeoTube("sV0ANail1PlaInHole7", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);    
    new TGeoTube("sV0ANail1PlaOuHole7", 0.0, 0.4, (fV0APlaAl)/2.); 

    /// Frame
    TGeoVolume *v0AFra7 = new TGeoVolumeAssembly("V0AFra7");
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = 0.;                              v0APts[1+8*i] = -fV0AR0-fV0AFraWd/2.;
    v0APts[2+8*i] = fV0AFraWd/2.;                    v0APts[3+8*i] = -fV0AR0-fV0AFraWd/2.;
    v0APts[4+8*i] = fV0AFraWd/2.;                    v0APts[5+8*i] = -fV0AR4+fV0AFraWd/2.;
    v0APts[6+8*i] = 0.;                              v0APts[7+8*i] = -fV0AR4+fV0AFraWd/2.;
    }
    TGeoArb8 *sV0AFraB17 = new TGeoArb8("sV0AFraB17",fV0ASciWd/2.,v0APts);
    TGeoVolume *v0AFraB17 = new TGeoVolume("V0AFraB17",sV0AFraB17,medV0AFra);
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = fV0AR0*cos45-fV0AFraWd;
    v0APts[1+8*i] =  -(fV0AR0-fV0AFraWd)*sin45;
    v0APts[2+8*i] = (fV0AR0-fV0AFraWd/2.)*cos45;
    v0APts[3+8*i] = -(fV0AR0-fV0AFraWd/2.)*sin45;
    v0APts[4+8*i] = (fV0AR4+fV0AFraWd/2.)*cos45/2.;
    v0APts[5+8*i] = -fV0AR4*sin45/2.;
    v0APts[6+8*i] = (fV0AR4+fV0AFraWd/4.)*cos45/2.;
    v0APts[7+8*i] = -(fV0AR4+fV0AFraWd/2.)*sin45/2.;
    }
    TGeoArb8 *sV0AFraB27 = new TGeoArb8("sV0AFraB27", fV0ASciWd/2., v0APts);
    TGeoVolume *v0AFraB27 = new TGeoVolume("V0AFraB27",sV0AFraB27,medV0AFra);
    v0AFraB17->SetLineColor(kV0AColorFra); v0AFraB27->SetLineColor(kV0AColorFra);
    v0AFra7->AddNode(v0AFraB17,1);
    v0AFra7->AddNode(v0AFraB27,1);  // Prefer 2 GeoObjects insted of 3 GeoMovements
    new TGeoTubeSeg( "sV0AFraR1b7", fV0AR0-fV0AFraWd/2.,
    	     fV0AR0+fV0AFraWd/2., fV0ASciWd/2., 270.0, 315.0);
    new TGeoTubeSeg( "sV0AFraR2b7", fV0AR1-fV0AFraWd/2.,
    	     fV0AR1+fV0AFraWd/2., fV0ASciWd/2., 270.0, 315.0);
    new TGeoTubeSeg( "sV0AFraR3b7", fV0AR2-fV0AFraWd/2.,
    	     fV0AR2+fV0AFraWd/2., fV0ASciWd/2., 270.0, 315.0);
    new TGeoTubeSeg( "sV0AFraR4b7", fV0AR3-fV0AFraWd/2.,
    	     fV0AR3+fV0AFraWd/2., fV0ASciWd/2., 270.0, 315.0);
    new TGeoTubeSeg( "sV0AFraR5b7", fV0AR4-fV0AFraWd/2.,
    	     fV0AR4+fV0AFraWd/2., fV0ASciWd/2., 270.0, 315.0);
    TGeoCompositeShape *sV0AFraR17 = new TGeoCompositeShape("sV0AFraR17","sV0AFraR1b7-sV0ACha7");
    TGeoCompositeShape *sV0AFraR27 = new TGeoCompositeShape("sV0AFraR27","sV0AFraR2b7-sV0ACha7");
    TGeoCompositeShape *sV0AFraR37 = new TGeoCompositeShape("sV0AFraR37","sV0AFraR3b7-sV0ACha7");
    TGeoCompositeShape *sV0AFraR47 = new TGeoCompositeShape("sV0AFraR47","sV0AFraR4b7-sV0ACha7");
    TGeoCompositeShape *sV0AFraR57 = new TGeoCompositeShape("sV0AFraR57","sV0AFraR5b7-sV0ACha7");
    TGeoVolume *v0AFraR17 = new TGeoVolume("V0AFraR17",sV0AFraR17,medV0AFra);
    TGeoVolume *v0AFraR27 = new TGeoVolume("V0AFraR27",sV0AFraR27,medV0AFra);
    TGeoVolume *v0AFraR37 = new TGeoVolume("V0AFraR37",sV0AFraR37,medV0AFra);
    TGeoVolume *v0AFraR47 = new TGeoVolume("V0AFraR47",sV0AFraR47,medV0AFra);
    TGeoVolume *v0AFraR57 = new TGeoVolume("V0AFraR57",sV0AFraR57,medV0AFra);
    v0AFraR17->SetLineColor(kV0AColorFra); v0AFraR27->SetLineColor(kV0AColorFra);
    v0AFraR37->SetLineColor(kV0AColorFra); v0AFraR47->SetLineColor(kV0AColorFra);
    v0AFraR57->SetLineColor(kV0AColorFra);
    v0AFra7->AddNode(v0AFraR17,1);
    v0AFra7->AddNode(v0AFraR27,1);
    v0AFra7->AddNode(v0AFraR37,1);
    v0AFra7->AddNode(v0AFraR47,1);
    v0AFra7->AddNode(v0AFraR57,1);
    v0ASec7->AddNode(v0AFra7,1);

    /// Sensitive scintilator
    TGeoVolume *v0ASci7 = new TGeoVolumeAssembly("V0ASci7");
    new TGeoTubeSeg( "sV0AR1b7", fV0AR0+fV0AFraWd/2.,
    	     fV0AR1-fV0AFraWd/2., fV0ASciWd/2., 270.0, 315.0);
    new TGeoTubeSeg( "sV0AR2b7", fV0AR1+fV0AFraWd/2.,
    	     fV0AR2-fV0AFraWd/2., fV0ASciWd/2., 270.0, 315.0);
    new TGeoTubeSeg( "sV0AR3b7", fV0AR2+fV0AFraWd/2.,
    	     fV0AR3-fV0AFraWd/2., fV0ASciWd/2., 270.0, 315.0);
    new TGeoTubeSeg( "sV0AR4b7", fV0AR3+fV0AFraWd/2.,
    	     fV0AR4-fV0AFraWd/2., fV0ASciWd/2., 270.0, 315.0);
    TGeoCompositeShape *sV0AR17 = new TGeoCompositeShape("sV0AR17","sV0AR1b7-sV0ACha7");
    TGeoCompositeShape *sV0AR27 = new TGeoCompositeShape("sV0AR27","sV0AR2b7-sV0ACha7");
    TGeoCompositeShape *sV0AR37 = new TGeoCompositeShape("sV0AR37","sV0AR3b7-sV0ACha7");
    TGeoCompositeShape *sV0AR47 = new TGeoCompositeShape("sV0AR47","sV0AR4b7-sV0ACha7");
    TGeoVolume *v0L17 = new TGeoVolume("V0L17",sV0AR17,medV0ASci);
    TGeoVolume *v0L27 = new TGeoVolume("V0L27",sV0AR27,medV0ASci);
    TGeoVolume *v0L37 = new TGeoVolume("V0L37",sV0AR37,medV0ASci);
    TGeoVolume *v0L47 = new TGeoVolume("V0L47",sV0AR47,medV0ASci);
    v0L17->SetLineColor(kV0AColorSci); v0L27->SetLineColor(kV0AColorSci);
    v0L37->SetLineColor(kV0AColorSci); v0L47->SetLineColor(kV0AColorSci);
    v0ASci7->AddNode(v0L17,1);
    v0ASci7->AddNode(v0L27,1);
    v0ASci7->AddNode(v0L37,1);
    v0ASci7->AddNode(v0L47,1);
    
    // Bunch of fibers
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = -10.0;	      v0APts[1+8*i] = 13.1;
    v0APts[2+8*i] = 10.0;            v0APts[3+8*i] = 13.1;
    v0APts[4+8*i] = 8.0;            v0APts[5+8*i] = -29.0;
    v0APts[6+8*i] = -12.0;	      v0APts[7+8*i] = -12.0;
    }   
    new TGeoArb8("sV0AFib17", 0.01, v0APts);      
    rot = new TGeoRotation("rot");
    rot->RotateX(-2.0); 
    rot->RotateY(0.0);
    rot->RotateZ(248.0);
    TGeoCombiTrans *fib17pos = new TGeoCombiTrans("fib17pos", 40.0 - 3.3, -50.0 + 1.5, 0.5, rot);
    fib17pos->RegisterYourself();
    TGeoCompositeShape *sV0AFib17Hole = new TGeoCompositeShape("sV0AFib17Hole", "sV0AFib17:fib17pos-sV0AFicR57");
    TGeoVolume *v0AFib17Hole = new TGeoVolume("V0AFib17Hole",sV0AFib17Hole,medV0AFib);
    v0AFib17Hole->SetLineColor(kV0AColorFib);
    new TGeoArb8("sV0AFib27", 0.01, v0APts);      
    rot = new TGeoRotation("rot");
    rot->RotateX(2.0); 
    rot->RotateY(0.0);
    rot->RotateZ(248.0);
    TGeoCombiTrans *fib27pos = new TGeoCombiTrans("fib27pos", 40.0 - 3.3, -50.0 + 1.5, -0.5, rot);
    fib27pos->RegisterYourself();
    TGeoCompositeShape *sV0AFib27Hole = new TGeoCompositeShape("sV0AFib27Hole", "sV0AFib27:fib27pos-sV0AFicR57");
    TGeoVolume *v0AFib27Hole = new TGeoVolume("V0AFib27Hole",sV0AFib27Hole,medV0AFib);
    v0AFib27Hole->SetLineColor(kV0AColorFib);
    TGeoVolume *v0AFib7 = new TGeoVolumeAssembly("V0AFib7");
    v0AFib7->AddNode(v0AFib17Hole,1); 
    v0AFib7->AddNode(v0AFib27Hole,1);
    v0ASec7->AddNode(v0AFib7,1);

    /// Non-sensitive scintilator
    new TGeoTubeSeg("sV0AR5S27", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShape7, 270.0, 315.0);
    TGeoCompositeShape *sV0AR57 = new TGeoCompositeShape("V0AR57","sV0AR5S27 - sV0ACha7");
    TGeoVolume *v0AR57 = new TGeoVolume("V0AR57",sV0AR57,medV0ASci);
    v0AR57->SetLineColor(kV0AColorSci);
    v0ASci7->AddNode(v0AR57,1);
    v0ASec7->AddNode(v0ASci7,1);

    /// Segment of octagon   
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = 0.;                 	 v0APts[1+8*i] = -fV0AR7*sin654;
    v0APts[2+8*i] = 0.;                          v0APts[3+8*i] = -(fV0AR7-fV0AOctH2)*sin654;
    v0APts[4+8*i] = (fV0AR7-fV0AOctH2)*cos654;   v0APts[5+8*i] = -(fV0AR7-fV0AOctH2)*sin654;
    v0APts[6+8*i] = fV0AR7*cos654;	         v0APts[7+8*i] = -fV0AR7*sin654;
    }
    TGeoArb8 *sV0AOct27 = new TGeoArb8("sV0AOct27", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoVolume *v0AOct27 = new TGeoVolume("V0AOct27", sV0AOct27,medV0ASup);
    v0AOct27->SetLineColor(kV0AColorOct);
    TGeoVolume *v0ASup7 = new TGeoVolumeAssembly("V0ASup7");
    v0ASup7->AddNode(v0AOct27,1);
    v0ASec7->AddNode(v0ASup7,1);

    /// Plates
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = 0.;       		v0APts[1+8*i] = -fV0AR7*sin654;
    v0APts[2+8*i] = 0.;        		v0APts[3+8*i] = -fV0AR0;
    v0APts[4+8*i] = fV0AR0*cos654;   	v0APts[5+8*i] = -fV0AR0*sin654;
    v0APts[6+8*i] = fV0AR7*cos654;  	v0APts[7+8*i] = -fV0AR7*sin654;
    }
    new TGeoArb8("sV0APlaIn7", (fV0APlaWd-2*fV0APlaAl)/2., v0APts);
    TGeoCompositeShape *sV0APlaInNailsHoles7 = new TGeoCompositeShape("sV0APlaInNailsHoles7","sV0APlaIn7-sV0ANail1PlaInHole7:pos17");
    TGeoVolume *v0APlaInNailsHoles7 = new TGeoVolume("V0APlaInNailsHoles7", sV0APlaInNailsHoles7, medV0APlaIn);
    new TGeoArb8("sV0APlaOu7", fV0APlaAl/2., v0APts);
    TGeoCompositeShape *sV0APlaOuNailsHoles7 = new TGeoCompositeShape("sV0APlaOuNailsHoles7","sV0APlaOu7-sV0ANail1PlaOuHole7:pos17"); 
    TGeoVolume *v0APlaOuNailsHoles7 = new TGeoVolume("V0APlaOuNailsHoles7", sV0APlaOuNailsHoles7, medV0APlaOu);
    v0APlaInNailsHoles7->SetLineColor(kV0AColorPlaIn); v0APlaOuNailsHoles7->SetLineColor(kV0AColorPlaOu);
    TGeoVolume *v0APla7 = new TGeoVolumeAssembly("V0APla7");
    v0APla7->AddNode(v0APlaInNailsHoles7,1);
    v0APla7->AddNode(v0APlaOuNailsHoles7,1,new TGeoTranslation(0,0,(fV0APlaWd-fV0APlaAl)/2.));
    v0APla7->AddNode(v0APlaOuNailsHoles7,2,new TGeoTranslation(0,0,-(fV0APlaWd-fV0APlaAl)/2.));
    v0ASec7->AddNode(v0APla7,1,new TGeoTranslation(0,0,(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    v0ASec7->AddNode(v0APla7,2,new TGeoTranslation(0,0,-(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    
    /// PMBox  
    TGeoVolume* v0APM7 = new TGeoVolumeAssembly("V0APM7");
    new TGeoBBox("sV0APMB17", fV0APMBWd/2., fV0APMBHt/2., fV0APMBTh/2.);
    new TGeoBBox("sV0APMB27", fV0APMBWd/2.-fV0APMBWdW, fV0APMBHt/2.-fV0APMBHtW, fV0APMBTh/2.-fV0APMBThW);
    TGeoCompositeShape *sV0APMB7 = new TGeoCompositeShape("sV0APMB7","sV0APMB17-sV0APMB27");
    TGeoVolume *v0APMB7 = new TGeoVolume("V0APMB7",sV0APMB7, medV0APMAlum);
    v0APMB7->SetLineColor(kV0AColorPMA);
    v0APM7->AddNode(v0APMB7,1);
    
    /// PMTubes 
    TGeoTube *sV0APMT17 = new TGeoTube("sV0APMT17", fV0APMTR1, fV0APMTR2, fV0APMTH/2.);
    TGeoVolume *v0APMT17 = new TGeoVolume("V0APMT17", sV0APMT17, medV0APMGlass);
    TGeoTube *sV0APMT27 = new TGeoTube("sV0APMT27", fV0APMTR3, fV0APMTR4, fV0APMTH/2.);
    TGeoVolume *v0APMT27 = new TGeoVolume("V0APMT27", sV0APMT27, medV0APMAlum);
    TGeoVolume *v0APMT7 = new TGeoVolumeAssembly("V0APMT7"); // pk si no choca con la 752 o con la 794
    TGeoTube *sV0APMTT7 = new TGeoTube("sV0APMTT7", 0., fV0APMTR4, fV0APMTB/2.);
    TGeoVolume *v0APMTT7 = new TGeoVolume("V0APMTT7", sV0APMTT7, medV0APMAlum);
    v0APMT7->SetLineColor(kV0AColorPMG);
    v0APMT27->SetLineColor(kV0AColorPMA);
    v0APMTT7->SetLineColor(kV0AColorPMA);
    rot = new TGeoRotation("rot", 90, 0, 180, 0, 90, 90);
    v0APMT7->AddNode(v0APMT17,1,rot);
    v0APMT7->AddNode(v0APMT27,1,rot);
    v0APMT7->AddNode(v0APMTT7,1,new TGeoCombiTrans(0,(fV0APMTH+fV0APMTB)/2.,0,rot));
    v0APM7->AddNode(v0APMT7, 1, new TGeoTranslation(-1.5*autoShiftSec1, 0, 0));
    v0APM7->AddNode(v0APMT7, 2, new TGeoTranslation(-0.5*autoShiftSec1, 0, 0));
    v0APM7->AddNode(v0APMT7, 3, new TGeoTranslation(+0.5*autoShiftSec1, 0, 0));
    v0APM7->AddNode(v0APMT7, 4, new TGeoTranslation(+1.5*autoShiftSec1, 0, 0));

    /// PM 
    rot = new TGeoRotation("rot");
    rot->RotateX(-90+30);
    rot->RotateY(0);
    rot->RotateZ(65+3);
    double shiftZ7 = fV0APMBHt/2. * cosAngPMBSec1
      -   ( fV0ASciWd + 2 * fV0AOctWd + 2 * fV0APlaWd )/2.   -   fV0APMBTh/2. * sinAngPMBSec1;
    double shiftR7 = fV0AR6  + fV0AR1 + fV0AOctWd + fV0APlaAl/3.;
    v0ASec7->AddNode(v0APM7,1, new TGeoCombiTrans( shiftR7*cos45+1.3, -shiftR7*sin45, shiftZ7, rot)); 
    
    // Aluminium nail
    TGeoTube *sV0ANail71 = new TGeoTube("sV0ANail71", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail71 = new TGeoVolume("V0ANail71", sV0ANail71, medV0APMAlum);
    v0ANail71->SetLineColor(kV0AColorPMA);// this is the color for aluminium
    v0ASec7->AddNode(v0ANail71,1,new TGeoTranslation(0.51,-42.9,0.0));
    TGeoTube *sV0ANail72 = new TGeoTube("sV0ANail72", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail72 = new TGeoVolume("V0ANail72", sV0ANail72, medV0APMAlum);
    v0ANail72->SetLineColor(kV0AColorPMA);
    v0ASec7->AddNode(v0ANail72,1,new TGeoTranslation(30.05,-30.79,0.0));  

    // Support
    TGeoBBox *sV0ASuppbr = new TGeoBBox("sV0ASuppbr", 2.0, 18.13, fV0ASciWd/2.);      
    TGeoVolume *v0ASuppbr = new TGeoVolume("V0ASuppbr", sV0ASuppbr, medV0ASup);
    v0ASuppbr->SetLineColor(kV0AColorOct);
    v0ASec7->AddNode(v0ASuppbr,1,new TGeoTranslation(2.0,-63.64,0.0));
    
    // Adding sector to v0LE volume 
    v0LE->AddNode(v0ASec7,1);
    

   //Definition of sector 8
   TGeoVolume *v0ASec8 = new TGeoVolumeAssembly("V0ASec8"); 
  
  /// For boolean sustraction
      double preShape8 = 0.2;
      for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[1+8*i] = -2*preShape8;
      v0APts[2+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[3+8*i] = 0.0;
      v0APts[4+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[5+8*i] = 0.0;
      v0APts[6+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[7+8*i] = -2*preShape8;
    }
    new TGeoArb8("sV0ACha18",fV0ASciWd/1.5,v0APts);
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = (fV0AR0-fV0AFraWd/2.)*sin45-preShape8;
      v0APts[1+8*i] = -(fV0AR0-fV0AFraWd/2.)*sin45;
      v0APts[2+8*i] = fV0AR0*sin45-fV0AFraWd/2.-preShape8;
      v0APts[3+8*i] = -(fV0AR0-fV0AFraWd)*sin45+preShape8;
      v0APts[4+8*i] = (fV0AR4+3*fV0AFraWd/2.)*sin45+preShape8;
      v0APts[5+8*i] = -fV0AR4*sin45-preShape8;
      v0APts[6+8*i] = (fV0AR4+fV0AFraWd)*sin45+preShape8; 
      v0APts[7+8*i] = -(fV0AR4+fV0AFraWd/2.)*sin45-2.*preShape8;
    }
    new TGeoArb8("sV0ACha28", fV0ASciWd/2.+2.*preShape8, v0APts);
    new TGeoCompositeShape("sV0ACha128","sV0ACha18+sV0ACha28");
    new TGeoTube("sV0ANail18Hole", 0.0, 0.4, 1.65);
    TGeoTranslation *pos18 = new TGeoTranslation("pos18",42.9,-.51,0.0);
    pos18->RegisterYourself();
    new TGeoTube("sV0ANail28Hole", 0.0, 0.4, 1.65);
    TGeoTranslation *pos28 = new TGeoTranslation("pos28",30.8,-30.04,0.0);
    pos28->RegisterYourself();    
    TGeoTranslation *pos38 = new TGeoTranslation("pos38",30.05,-30.79,0.0);  
    pos38->RegisterYourself();
    new TGeoCompositeShape("sV0ANailsHoles8","sV0ANail18Hole:pos18+sV0ANail28Hole:pos28");
    new TGeoCompositeShape("sV0ACha8","sV0ACha128+sV0ANailsHoles8");
    new TGeoTubeSeg("sV0AFicR58", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2. +2*preShape8, 315.0, 360.0);
    new TGeoTube("sV0ANail1PlaInHole8", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoTube("sV0ANail2PlaInHole8", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoTube("sV0ANail3PlaInHole8", 0.0, 0.4, (fV0APlaWd-2*fV0APlaAl)/2.);
    new TGeoCompositeShape("sV0ANailsPlaInHoles8","sV0ANail1PlaInHole8:pos18+sV0ANail2PlaInHole8:pos28+sV0ANail3PlaInHole8:pos38");
    new TGeoTube("sV0ANail1PlaOuHole8", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoTube("sV0ANail2PlaOuHole8", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoTube("sV0ANail3PlaOuHole8", 0.0, 0.4, (fV0APlaAl)/2.);
    new TGeoCompositeShape("sV0ANailsPlaOuHoles8","sV0ANail1PlaOuHole8:pos18+sV0ANail2PlaOuHole8:pos28+sV0ANail3PlaOuHole8:pos38");
    
    /// Frame 
    TGeoVolume *v0AFra8 = new TGeoVolumeAssembly("V0AFra8"); 
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[1+8*i] = 0.0;
      v0APts[2+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[3+8*i] = 0.0;
      v0APts[4+8*i] = fV0AR4+fV0AFraWd/2.;  v0APts[5+8*i] = -fV0AFraWd/8.;
      v0APts[6+8*i] = fV0AR0-fV0AFraWd/2.;  v0APts[7+8*i] = -fV0AFraWd/8.;
    }    
    TGeoArb8 *sV0AFraB18 = new TGeoArb8("sV0AFraB18",fV0ASciWd/2.,v0APts);
    TGeoVolume *v0AFraB18 = new TGeoVolume("V0AFraB18",sV0AFraB18,medV0AFra);  
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = (fV0AR0-fV0AFraWd/4.)*sin45;
      v0APts[1+8*i] = -(fV0AR0-fV0AFraWd/2.)*sin45; 
      v0APts[2+8*i] = fV0AR0*sin45-fV0AFraWd/2.; 
      v0APts[3+8*i] = -(fV0AR0-fV0AFraWd)*sin45;
      v0APts[4+8*i] = (fV0AR4+3*fV0AFraWd/2.)*sin45/2.;
      v0APts[5+8*i] = -fV0AR4*sin45/2.;
      v0APts[6+8*i] = (fV0AR4+fV0AFraWd)*sin45/2.;
      v0APts[7+8*i] = -(fV0AR4+fV0AFraWd/2.)*sin45/2.;
    }
    TGeoArb8 *sV0AFraB28 = new TGeoArb8("sV0AFraB28", fV0ASciWd/2., v0APts);
    TGeoVolume *v0AFraB28 = new TGeoVolume("V0AFraB28",sV0AFraB28,medV0AFra);
    v0AFraB18->SetLineColor(kV0AColorFra); v0AFraB28->SetLineColor(kV0AColorFra);
    v0AFra8->AddNode(v0AFraB18,1);
    v0AFra8->AddNode(v0AFraB28,1);  // Prefer 2 GeoObjects insted of 3 GeoMovements
    new TGeoTubeSeg( "sV0AFraR1b8", fV0AR0-fV0AFraWd/2.,
		     fV0AR0+fV0AFraWd/2., fV0ASciWd/2., 315.0, 360.0);
    new TGeoTubeSeg( "sV0AFraR2b8", fV0AR1-fV0AFraWd/2.,
		     fV0AR1+fV0AFraWd/2., fV0ASciWd/2., 315.0, 360.0);
    new TGeoTubeSeg( "sV0AFraR3b8", fV0AR2-fV0AFraWd/2.,
		     fV0AR2+fV0AFraWd/2., fV0ASciWd/2., 315.0, 360.0);
    new TGeoTubeSeg( "sV0AFraR4b8", fV0AR3-fV0AFraWd/2.,
		     fV0AR3+fV0AFraWd/2., fV0ASciWd/2., 315.0, 360.0);
    new TGeoTubeSeg( "sV0AFraR5b8", fV0AR4-fV0AFraWd/2.,
		     fV0AR4+fV0AFraWd/2., fV0ASciWd/2., 315.0, 360.0);
    TGeoCompositeShape *sV0AFraR18 = new TGeoCompositeShape("sV0AFraR18","sV0AFraR1b8-sV0ACha8");
    TGeoCompositeShape *sV0AFraR28 = new TGeoCompositeShape("sV0AFraR28","sV0AFraR2b8-sV0ACha8");
    TGeoCompositeShape *sV0AFraR38 = new TGeoCompositeShape("sV0AFraR38","sV0AFraR3b8-sV0ACha8");
    TGeoCompositeShape *sV0AFraR48 = new TGeoCompositeShape("sV0AFraR48","sV0AFraR4b8-sV0ACha8");
    TGeoCompositeShape *sV0AFraR58 = new TGeoCompositeShape("sV0AFraR58","sV0AFraR5b8-sV0ACha8");
    TGeoVolume *v0AFraR18 = new TGeoVolume("V0AFraR18",sV0AFraR18,medV0AFra);
    TGeoVolume *v0AFraR28 = new TGeoVolume("V0AFraR28",sV0AFraR28,medV0AFra);
    TGeoVolume *v0AFraR38 = new TGeoVolume("V0AFraR38",sV0AFraR38,medV0AFra);
    TGeoVolume *v0AFraR48 = new TGeoVolume("V0AFraR48",sV0AFraR48,medV0AFra);
    TGeoVolume *v0AFraR58 = new TGeoVolume("V0AFraR58",sV0AFraR58,medV0AFra);
    v0AFraR18->SetLineColor(kV0AColorFra); v0AFraR28->SetLineColor(kV0AColorFra);
    v0AFraR38->SetLineColor(kV0AColorFra); v0AFraR48->SetLineColor(kV0AColorFra);
    v0AFraR58->SetLineColor(kV0AColorFra);
    v0AFra8->AddNode(v0AFraR18,1);
    v0AFra8->AddNode(v0AFraR28,1);
    v0AFra8->AddNode(v0AFraR38,1);
    v0AFra8->AddNode(v0AFraR48,1);
    v0AFra8->AddNode(v0AFraR58,1);
    v0ASec8->AddNode(v0AFra8,1);

    /// Sensitive scintilator
    TGeoVolume *v0ASci8 = new TGeoVolumeAssembly("V0ASci8");
    new TGeoTubeSeg( "sV0AR1b8", fV0AR0+fV0AFraWd/2.,
		     fV0AR1-fV0AFraWd/2., fV0ASciWd/2., 315.0, 360.0);
    new TGeoTubeSeg( "sV0AR2b8", fV0AR1+fV0AFraWd/2.,
		     fV0AR2-fV0AFraWd/2., fV0ASciWd/2., 315.0, 360.0);
    new TGeoTubeSeg( "sV0AR3b8", fV0AR2+fV0AFraWd/2.,
		     fV0AR3-fV0AFraWd/2., fV0ASciWd/2., 315.0, 360.0);
    new TGeoTubeSeg( "sV0AR4b8", fV0AR3+fV0AFraWd/2.,
		     fV0AR4-fV0AFraWd/2., fV0ASciWd/2., 315.0, 360.0);
    TGeoCompositeShape *sV0AR18 = new TGeoCompositeShape("sV0AR18","sV0AR1b8-sV0ACha8");
    TGeoCompositeShape *sV0AR28 = new TGeoCompositeShape("sV0AR28","sV0AR2b8-sV0ACha8");
    TGeoCompositeShape *sV0AR38 = new TGeoCompositeShape("sV0AR38","sV0AR3b8-sV0ACha8");
    TGeoCompositeShape *sV0AR48 = new TGeoCompositeShape("sV0AR48","sV0AR4b8-sV0ACha8");
    TGeoVolume *v0L18 = new TGeoVolume("V0L18",sV0AR18,medV0ASci);
    TGeoVolume *v0L28 = new TGeoVolume("V0L28",sV0AR28,medV0ASci);
    TGeoVolume *v0L38 = new TGeoVolume("V0L38",sV0AR38,medV0ASci);
    TGeoVolume *v0L48 = new TGeoVolume("V0L48",sV0AR48,medV0ASci);
    v0L18->SetLineColor(kV0AColorSci); v0L28->SetLineColor(kV0AColorSci);
    v0L38->SetLineColor(kV0AColorSci); v0L48->SetLineColor(kV0AColorSci);
    v0ASci8->AddNode(v0L18,1);
    v0ASci8->AddNode(v0L28,1);
    v0ASci8->AddNode(v0L38,1);
    v0ASci8->AddNode(v0L48,1); 

    /// Segment of octagon   
    for (int i=0;i<2;i++) {
    v0APts[0+8*i] = fV0AR6;	                 v0APts[1+8*i] = 0.;
    v0APts[2+8*i] = fV0AR7*cos654;               v0APts[3+8*i] = -fV0AR7*sin654;
    v0APts[4+8*i] = (fV0AR7-fV0AOctH2)*cos654;   v0APts[5+8*i] = -(fV0AR7-fV0AOctH2)*sin654;
    v0APts[6+8*i] = fV0AR6-fV0AOctH2;		 v0APts[7+8*i] = 0.;
    }
    new TGeoArb8("sV0AOct28", (fV0ASciWd+2*fV0AOctWd)/2., v0APts);
    TGeoCompositeShape *sV0AOct2FEEB8 = new TGeoCompositeShape("sV0AOct2FEEB8","sV0AOct28-sV0AFicFEEBSec1:posFicFEEBSec1-sV0AFicOct2Sec1:posFicOct2Sec1");
    TGeoVolume *v0AOct28 = new TGeoVolume("V0AOct28", sV0AOct2FEEB8,medV0ASup);
    v0AOct28->SetLineColor(kV0AColorOct);
    TGeoVolume *v0ASup8 = new TGeoVolumeAssembly("V0ASup8");
    v0ASup8->AddNode(v0AOct28,1);
    v0ASec8->AddNode(v0ASup8,1);

    //Bunch of fibers
    v0APts[ 0] = v0APts[ 2] = -14.0;
    v0APts[ 1] = v0APts[ 7] = (fV0ASciWd+fV0AOctWd)/2.-0.01;
    v0APts[ 3] = v0APts[ 5] = (fV0ASciWd+fV0AOctWd)/2.+0.01;
    v0APts[ 4] = v0APts[ 6] = +14.0;
    v0APts[ 8] = v0APts[10] = -10.0;
    v0APts[ 9] = v0APts[15] = 0.;
    v0APts[11] = v0APts[13] = 0.25;
    v0APts[12] = v0APts[14] = +10.0;
    new TGeoArb8("sV0AFib18", 11.8, v0APts); 
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateZ(-90-22.5);
    TGeoCombiTrans *fib18pos = new TGeoCombiTrans("fib18pos", (fV0AR6-fV0AOctH2+fV0AR5)*cos225/2. - 3.3, -(fV0AR6-fV0AOctH2+fV0AR5)*sin225/2. + 1.5, 0, rot);
    fib18pos->RegisterYourself();
    TGeoCompositeShape *sV0AFib18Hole = new TGeoCompositeShape("sV0AFib18Hole", "sV0AFib18:fib18pos-sV0AFicR58");
    TGeoVolume *v0AFib18Hole = new TGeoVolume("V0AFib18",sV0AFib18Hole,medV0AFib);
    v0AFib18Hole->SetLineColor(kV0AColorFib);
    new TGeoArb8("sV0AFib28", 11.8, v0APts);
    rot = new TGeoRotation("rot");
    rot->RotateX(-90);
    rot->RotateY(180);
    rot->RotateZ(-90-22.5);
    TGeoCombiTrans *fib28pos = new TGeoCombiTrans("fib28pos", (fV0AR6-fV0AOctH2+fV0AR5)*cos225/2. - 3.3, -(fV0AR6-fV0AOctH2+fV0AR5)*sin225/2. + 1.5, 0, rot);
    fib28pos->RegisterYourself();
    TGeoCompositeShape *sV0AFib28Hole = new TGeoCompositeShape("sV0AFib28Hole", "sV0AFib28:fib28pos-sV0AFicR58");
    TGeoVolume *v0AFib28Hole = new TGeoVolume("V0AFib28Hole",sV0AFib28Hole,medV0AFib);
    v0AFib28Hole->SetLineColor(kV0AColorFib);
    TGeoVolume *v0AFib8 = new TGeoVolumeAssembly("V0AFib8");    
    v0AFib8->AddNode(v0AFib18Hole,1);
    v0AFib8->AddNode(v0AFib28Hole,1);
    v0ASec8->AddNode(v0AFib8,1);
             
    /// Non-sensitive scintilator   
    new TGeoTubeSeg("sV0AR5S28", fV0AR4+fV0AFraWd/2., fV0AR4 + fV0AR0, fV0ASciWd/2.+2*preShape8, 315.0, 360.0);
    TGeoCompositeShape *sV0AR58 = new TGeoCompositeShape("V0AR58","sV0AR5S28 - sV0ACha8");
    TGeoVolume *v0AR58 = new TGeoVolume("V0AR58",sV0AR58,medV0ASci);
    v0AR58->SetLineColor(kV0AColorSci);
    v0ASci8->AddNode(v0AR58,1);
    v0ASec8->AddNode(v0ASci8,1); 
    
    /// Plates
    for (int i=0;i<2;i++) {
      v0APts[0+8*i] = fV0AR0;			v0APts[1+8*i] = 0.;
      v0APts[2+8*i] = fV0AR6;   		v0APts[3+8*i] = 0.;
      v0APts[4+8*i] = fV0AR7*cos654;    	v0APts[5+8*i] = -fV0AR7*sin654;
      v0APts[6+8*i] = fV0AR0*cos654;  		v0APts[7+8*i] = -fV0AR0*sin654;
    }
    new TGeoArb8("sV0APlaIn8", (fV0APlaWd-2*fV0APlaAl)/2., v0APts);
    TGeoCompositeShape *sV0APlaInNailsHoles8 = new TGeoCompositeShape("sV0APlaInNailsHoles8","sV0APlaIn8-sV0ANailsPlaInHoles8-sV0AFicFEEBSec1:posFicFEEBSec1");
    TGeoVolume *v0APlaInNailsHoles8 = new TGeoVolume("V0APlaInNailsHoles8", sV0APlaInNailsHoles8, medV0APlaIn);
    new TGeoArb8("sV0APlaOu8", fV0APlaAl/2., v0APts);
    TGeoCompositeShape *sV0APlaOuNailsHoles8 = new TGeoCompositeShape("sV0APlaOuNailsHoles8","sV0APlaOu8-sV0ANailsPlaOuHoles8-sV0AFicFEEBSec1:posFicFEEBSec1"); 
    TGeoVolume *v0APlaOuNailsHoles8 = new TGeoVolume("V0APlaOuNailsHoles8", sV0APlaOuNailsHoles8, medV0APlaOu);
    v0APlaInNailsHoles8->SetLineColor(kV0AColorPlaIn); v0APlaOuNailsHoles8->SetLineColor(kV0AColorPlaOu);
    TGeoVolume *v0APla8 = new TGeoVolumeAssembly("V0APla8");
    v0APla8->AddNode(v0APlaInNailsHoles8,1);
    v0APla8->AddNode(v0APlaOuNailsHoles8,1,new TGeoTranslation(0,0,(fV0APlaWd-fV0APlaAl)/2.));
    v0APla8->AddNode(v0APlaOuNailsHoles8,2,new TGeoTranslation(0,0,-(fV0APlaWd-fV0APlaAl)/2.));
    v0ASec8->AddNode(v0APla8,1,new TGeoTranslation(0,0,(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));
    v0ASec8->AddNode(v0APla8,2,new TGeoTranslation(0,0,-(fV0ASciWd+2*fV0AOctWd+fV0APlaWd)/2.));

    /// PMBox 
    TGeoVolume* v0APM8 = new TGeoVolumeAssembly("V0APM1");
    new TGeoBBox("sV0APMB18", fV0APMBWd/2., fV0APMBHt/2., fV0APMBTh/2.);
    new TGeoBBox("sV0APMB28", fV0APMBWd/2.-fV0APMBWdW, fV0APMBHt/2.-fV0APMBHtW, fV0APMBTh/2.-fV0APMBThW);
    TGeoCompositeShape *sV0APMB8 = new TGeoCompositeShape("sV0APMB8","sV0APMB18-sV0APMB28");
    TGeoVolume *v0APMB8 = new TGeoVolume("V0APMB8",sV0APMB8, medV0APMAlum);
    v0APMB8->SetLineColor(kV0AColorPMA);
    v0APM8->AddNode(v0APMB8,1);

    /// PMTubes 
    TGeoTube *sV0APMT18 = new TGeoTube("sV0APMT18", fV0APMTR1, fV0APMTR2, fV0APMTH/2.);
    TGeoVolume *v0APMT18 = new TGeoVolume("V0APMT18", sV0APMT18, medV0APMGlass);
    TGeoTube *sV0APMT28 = new TGeoTube("sV0APMT28", fV0APMTR3, fV0APMTR4, fV0APMTH/2.);
    TGeoVolume *v0APMT28 = new TGeoVolume("V0APMT28", sV0APMT28, medV0APMAlum);
    TGeoVolume *v0APMT8 = new TGeoVolumeAssembly("V0APMT8");
    TGeoTube *sV0APMTT8 = new TGeoTube("sV0APMTT8", 0., fV0APMTR4, fV0APMTB/2.);
    TGeoVolume *v0APMTT8 = new TGeoVolume("V0APMTT8", sV0APMTT8, medV0APMAlum);
    v0APMT8->SetLineColor(kV0AColorPMG);
    v0APMT28->SetLineColor(kV0AColorPMA);
    v0APMTT8->SetLineColor(kV0AColorPMA);
    rot = new TGeoRotation("rot", 90, 0, 180, 0, 90, 90);
    v0APMT8->AddNode(v0APMT18,1,rot);
    v0APMT8->AddNode(v0APMT28,1,rot);
    v0APMT8->AddNode(v0APMTT8,1,new TGeoCombiTrans(0,(fV0APMTH+fV0APMTB)/2.,0,rot));
    double autoShift8 = (fV0APMBWd-2*fV0APMBWdW)/4.;
    v0APM8->AddNode(v0APMT8, 1, new TGeoTranslation(-1.5*autoShift8, 0, 0));
    v0APM8->AddNode(v0APMT8, 2, new TGeoTranslation(-0.5*autoShift8, 0, 0));
    v0APM8->AddNode(v0APMT8, 3, new TGeoTranslation(+0.5*autoShift8, 0, 0));
    v0APM8->AddNode(v0APMT8, 4, new TGeoTranslation(+1.5*autoShift8, 0, 0));

    /// PM 
    rot = new TGeoRotation("rot");
    rot->RotateX(-90+30);
    rot->RotateY(0);
    rot->RotateZ(65+3);
    double shiftZ8 = fV0APMBHt/2. * cosAngPMBSec1
      -   ( fV0ASciWd + 2 * fV0AOctWd + 2 * fV0APlaWd )/2.   -   fV0APMBTh/2. * sinAngPMBSec1;
    double shiftR8 = fV0AR6  +  fV0AOctH2 + fV0APlaAl;
    v0ASec8->AddNode(v0APM8,1, new TGeoCombiTrans( shiftR8*cos225 + 1.3, -shiftR8*sin225, shiftZ8, rot));

    // Aluminium nails
    TGeoTube *sV0ANail81 = new TGeoTube("sV0ANail81", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail81 = new TGeoVolume("V0ANail81", sV0ANail81, medV0APMAlum);
    v0ANail81->SetLineColor(kV0AColorPMA);// this is the color for aluminium
    v0ASec8->AddNode(v0ANail81,1,new TGeoTranslation(42.9,-.51,0.0));
    TGeoTube *sV0ANail82 = new TGeoTube("sV0ANail82", 0.0, 0.4, 5.09/2.);
    TGeoVolume *v0ANail82 = new TGeoVolume("V0ANail82", sV0ANail82, medV0APMAlum); 
    v0ANail82->SetLineColor(kV0AColorPMA);
    v0ASec8->AddNode(v0ANail82,1,new TGeoTranslation(30.8,-30.04,0.0));  
      
    // Adding sector to v0LE volume 
    v0LE->AddNode(v0ASec8, 1);    
      
    // Adding detectors to top volume
    TGeoVolume *vZERO = new TGeoVolumeAssembly("VZERO");
    vZERO->AddNode(v0RI,1,new TGeoTranslation(0, 0, -zdet));
    vZERO->AddNode(v0LE,1,new TGeoTranslation(0, 0, +329.0));
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
  Int_t     fieldType       = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();     // Field type 
  Double_t  maxField        = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();       // Field max.
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

// Parameters for V0AFiberGlass: Material for mother and daughter boards
   as[0] = 1.00794;	as[1] = 12.011;       as[2] = 16.0;	     as[3] = 28.09;
   zs[0] = 1.;		zs[1] = 6.;           zs[2] = 8.;	     zs[3] = 14.;
   ws[0] = 736.0;	ws[1] = 462.0;        ws[2] = 292.0;	     ws[3] = 68.0;
   density      = 1.9;
   id           = 12;
   AliMixture(id, "V0AFibGlass", as, zs, density, -4, ws);
   AliMedium(id, "V0AFibGlass", id, 1, fieldType, maxField, maxBending, maxStepSize,
	     maxEnergyLoss, precision, minStepSize); 
   
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
  static    Int_t   nPhotonsInStep = 0;
  static    Int_t   nPhotons = 0; 
  static    Int_t   numStep = 0;
  Int_t     ringNumber;
  Float_t   destep, step;
  numStep += 1; 

  //   We keep only charged tracks : 
  if ( !fMC->TrackCharge() || !fMC->IsTrackAlive() ) return;

  vol[0]    = fMC->CurrentVolOffID(1, vol[1]);
  vol[2]    = fMC->CurrentVolID(copy);
  vol[3]    = copy;
  static Int_t idV0R1 = fMC->VolId("V0R1");

  static Int_t idV0L11 = fMC->VolId("V0L1Sec1");
  static Int_t idV0L12 = fMC->VolId("V0L1Sec2");
  static Int_t idV0L13 = fMC->VolId("V0L1Sec3");
  static Int_t idV0L14 = fMC->VolId("V0L1Sec4");
  static Int_t idV0L15 = fMC->VolId("V0L15");
  static Int_t idV0L16 = fMC->VolId("V0L16");
  static Int_t idV0L17 = fMC->VolId("V0L17");
  static Int_t idV0L18 = fMC->VolId("V0L18");
  static Int_t idV0R2 = fMC->VolId("V0R2");

  static Int_t idV0L21 = fMC->VolId("V0L2Sec1");
  static Int_t idV0L22 = fMC->VolId("V0L2Sec2");
  static Int_t idV0L23 = fMC->VolId("V0L2Sec3");
  static Int_t idV0L24 = fMC->VolId("V0L2Sec4");
  static Int_t idV0L25 = fMC->VolId("V0L25");
  static Int_t idV0L26 = fMC->VolId("V0L26");
  static Int_t idV0L27 = fMC->VolId("V0L27");
  static Int_t idV0L28 = fMC->VolId("V0L28");
  static Int_t idV0R3 = fMC->VolId("V0R3");

  static Int_t idV0L31 = fMC->VolId("V0L3Sec1");
  static Int_t idV0L32 = fMC->VolId("V0L3Sec2");
  static Int_t idV0L33 = fMC->VolId("V0L3Sec3");
  static Int_t idV0L34 = fMC->VolId("V0L3Sec4");
  static Int_t idV0L35 = fMC->VolId("V0L35");
  static Int_t idV0L36 = fMC->VolId("V0L36");
  static Int_t idV0L37 = fMC->VolId("V0L37");
  static Int_t idV0L38 = fMC->VolId("V0L38");
  static Int_t idV0R4 = fMC->VolId("V0R4");

  static Int_t idV0L41 = fMC->VolId("V0L4Sec1");
  static Int_t idV0L42 = fMC->VolId("V0L4Sec2");
  static Int_t idV0L43 = fMC->VolId("V0L4Sec3");
  static Int_t idV0L44 = fMC->VolId("V0L4Sec4");
  static Int_t idV0L45 = fMC->VolId("V0L45");
  static Int_t idV0L46 = fMC->VolId("V0L46");
  static Int_t idV0L47 = fMC->VolId("V0L47");
  static Int_t idV0L48 = fMC->VolId("V0L48");
  static Int_t idV0R5 = fMC->VolId("V0R5");
  static Int_t idV0R6 = fMC->VolId("V0R6");
  bool   hitOnV0C = true;
  double lightYield;
  double lightAttenuation;
  double nMeters; 
  double fibToPhot;
  if      ( fMC->CurrentVolID(copy) == idV0R1   ||
        fMC->CurrentVolID(copy) == idV0L11  ||
        fMC->CurrentVolID(copy) == idV0L12  ||
        fMC->CurrentVolID(copy) == idV0L13  ||
        fMC->CurrentVolID(copy) == idV0L14  ||
        fMC->CurrentVolID(copy) == idV0L15  ||
        fMC->CurrentVolID(copy) == idV0L16  ||
        fMC->CurrentVolID(copy) == idV0L17  ||
        fMC->CurrentVolID(copy) == idV0L18
      )
      ringNumber = 1;
  
  else if ( fMC->CurrentVolID(copy) == idV0R2  ||
        fMC->CurrentVolID(copy) == idV0L21 ||
        fMC->CurrentVolID(copy) == idV0L22 ||
        fMC->CurrentVolID(copy) == idV0L23 ||
        fMC->CurrentVolID(copy) == idV0L24 ||
        fMC->CurrentVolID(copy) == idV0L25 ||
        fMC->CurrentVolID(copy) == idV0L26 ||
        fMC->CurrentVolID(copy) == idV0L27 ||
        fMC->CurrentVolID(copy) == idV0L28
      )
      ringNumber = 2; 
  
  else if ( fMC->CurrentVolID(copy) == idV0R3  ||
        fMC->CurrentVolID(copy) == idV0R4  ||
        fMC->CurrentVolID(copy) == idV0L31 ||
        fMC->CurrentVolID(copy) == idV0L32 ||
        fMC->CurrentVolID(copy) == idV0L33 ||
        fMC->CurrentVolID(copy) == idV0L34 ||
        fMC->CurrentVolID(copy) == idV0L35 ||
        fMC->CurrentVolID(copy) == idV0L36 ||
        fMC->CurrentVolID(copy) == idV0L37 ||
        fMC->CurrentVolID(copy) == idV0L38
      ) 
      ringNumber = 3;
  else if ( fMC->CurrentVolID(copy) == idV0R5  ||
        fMC->CurrentVolID(copy) == idV0R6  ||
        fMC->CurrentVolID(copy) == idV0L41 ||
        fMC->CurrentVolID(copy) == idV0L42 ||
        fMC->CurrentVolID(copy) == idV0L43 ||
        fMC->CurrentVolID(copy) == idV0L44 ||
        fMC->CurrentVolID(copy) == idV0L45 ||
        fMC->CurrentVolID(copy) == idV0L46 ||
        fMC->CurrentVolID(copy) == idV0L47 ||
        fMC->CurrentVolID(copy) == idV0L48
      ) ringNumber = 4;	       
  
  else ringNumber = 0;
  
  if  (ringNumber) {
      if (
      fMC->CurrentVolID(copy) == idV0L11 ||
      fMC->CurrentVolID(copy) == idV0L12 ||
      fMC->CurrentVolID(copy) == idV0L13 ||
      fMC->CurrentVolID(copy) == idV0L14 ||
      fMC->CurrentVolID(copy) == idV0L15 ||
      fMC->CurrentVolID(copy) == idV0L16 ||
      fMC->CurrentVolID(copy) == idV0L17 ||
      fMC->CurrentVolID(copy) == idV0L18 ||
      fMC->CurrentVolID(copy) == idV0L21 ||
      fMC->CurrentVolID(copy) == idV0L22 ||
      fMC->CurrentVolID(copy) == idV0L23 ||
      fMC->CurrentVolID(copy) == idV0L24 ||
      fMC->CurrentVolID(copy) == idV0L25 ||
      fMC->CurrentVolID(copy) == idV0L26 ||
      fMC->CurrentVolID(copy) == idV0L27 ||
      fMC->CurrentVolID(copy) == idV0L28 ||
      fMC->CurrentVolID(copy) == idV0L31 ||
      fMC->CurrentVolID(copy) == idV0L32 ||
      fMC->CurrentVolID(copy) == idV0L33 ||
      fMC->CurrentVolID(copy) == idV0L34 ||
      fMC->CurrentVolID(copy) == idV0L35 ||
      fMC->CurrentVolID(copy) == idV0L36 ||
      fMC->CurrentVolID(copy) == idV0L37 ||
      fMC->CurrentVolID(copy) == idV0L38 ||
      fMC->CurrentVolID(copy) == idV0L41 ||
      fMC->CurrentVolID(copy) == idV0L42 ||
      fMC->CurrentVolID(copy) == idV0L43 ||
      fMC->CurrentVolID(copy) == idV0L44 ||
      fMC->CurrentVolID(copy) == idV0L45 ||
      fMC->CurrentVolID(copy) == idV0L46 ||
      fMC->CurrentVolID(copy) == idV0L47 ||
      fMC->CurrentVolID(copy) == idV0L48
	  )
	  hitOnV0C = false;

    destep = fMC->Edep();
    step   = fMC->TrackStep();
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
    if ( fMC->IsTrackEntering() ) {
      nPhotons  =  nPhotonsInStep;
      fMC->TrackPosition(fTrackPosition);
      fMC->TrackMomentum(fTrackMomentum);
      Float_t pt  = TMath::Sqrt( fTrackMomentum.Px() * fTrackMomentum.Px()
				 + fTrackMomentum.Py() * fTrackMomentum.Py() );
      TParticle *par = gAlice->GetMCApp()->Particle(gAlice->GetMCApp()->GetCurrentTrackNumber());
      hits[0]  = fTrackPosition.X();
      hits[1]  = fTrackPosition.Y();
      hits[2]  = fTrackPosition.Z();	 	 
      hits[3]  = Float_t (fMC->TrackPid());
      hits[4]  = fMC->TrackTime();
      hits[5]  = fMC->TrackCharge();
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
    if( fMC->IsTrackExiting() || fMC->IsTrackStop() || fMC->IsTrackDisappeared()){
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
    if( fMC->IsTrackEntering() || fMC->IsTrackExiting() ) {
      AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kVZERO);
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
void AliVZEROv7::MakeBranch(Option_t *option)
{
// Creates new branches in the current Root Tree
    
  TString branchname(Form("%s",GetName()));
  AliDebug(2,Form("fBufferSize = %d",fBufferSize));
  const char *cH = strstr(option,"H");
  if (fHits   && fLoader->TreeH() && cH) {
    fLoader->TreeH()->Branch(branchname.Data(),&fHits, fBufferSize);
    AliDebug(2,Form("Making Branch %s for hits",branchname.Data()));
  }     
  const char *cD = strstr(option,"D");
  if (fDigits   && fLoader->TreeD() && cD) {
    fLoader->TreeD()->Branch(branchname.Data(),&fDigits, fBufferSize);
    AliDebug(2,Form("Making Branch %s for digits",branchname.Data()));
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
      if ( TVirtualMC::GetMC()->CurrentVolID(vol[1]) == TVirtualMC::GetMC()->VolId("V0R3") || TVirtualMC::GetMC()->CurrentVolID(vol[1])
	   == TVirtualMC::GetMC()->VolId("V0R5") )  index = (index*2-14)+((ringNumber-2)*16);
      if ( TVirtualMC::GetMC()->CurrentVolID(vol[1]) == TVirtualMC::GetMC()->VolId("V0R4") || TVirtualMC::GetMC()->CurrentVolID(vol[1])
	   == TVirtualMC::GetMC()->VolId("V0R6") )  index = (index*2-13)+((ringNumber-2)*16);
    }
    fCellId   = index;           
  } else if (hits[2] > 0.0) {
    //    cout << " vol[0] = " << vol[0] << " copy : " << vol[1] 
    //	 << " called " << TVirtualMC::GetMC()->VolName(vol[0]) << endl;
    // cout << " vol[2] = " << vol[2] << " copy : " << vol[3] 
    //	 << " called " << TVirtualMC::GetMC()->VolName(vol[2]) << endl;
    // upper half

      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L1Sec1")) fCellId =  47 + 1;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L1Sec2")) fCellId =  47 + 2;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L1Sec3")) fCellId =  47 + 3;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L1Sec4")) fCellId =  47 + 4;

      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L2Sec1")) fCellId =  47 +  9;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L2Sec2")) fCellId =  47 + 10;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L2Sec3")) fCellId =  47 + 11;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L2Sec4")) fCellId =  47 + 12;

      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L3Sec1")) fCellId =  47 + 17;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L3Sec2")) fCellId =  47 + 18;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L3Sec3")) fCellId =  47 + 19;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L3Sec4")) fCellId =  47 + 20;

      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L4Sec1")) fCellId =  47 + 25;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L4Sec2")) fCellId =  47 + 26;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L4Sec3")) fCellId =  47 + 27;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L4Sec4")) fCellId =  47 + 28;

    // lower half 
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L15")) fCellId = 48+4;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L16")) fCellId = 48+5;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L17")) fCellId = 48+6;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L18")) fCellId = 48+7;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L25")) fCellId = 8+48+4;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L26")) fCellId = 8+48+5;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L27")) fCellId = 8+48+6;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L28")) fCellId = 8+48+7;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L35")) fCellId = 16+48+4;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L36")) fCellId = 16+48+5;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L37")) fCellId = 16+48+6;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L38")) fCellId = 16+48+7;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L45")) fCellId = 24+48+4;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L46")) fCellId = 24+48+5;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L47")) fCellId = 24+48+6;
      if (TVirtualMC::GetMC()->CurrentVolID(vol[2]) == TVirtualMC::GetMC()->VolId("V0L48")) fCellId = 24+48+7;
  }

  return fCellId;
}

