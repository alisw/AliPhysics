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

//-------------------------------------------------------------------------
//  Beam pipe class
//  This version uses TGeo
//  Author: A.Morsch
//-------------------------------------------------------------------------

#include <Riostream.h>

#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoTorus.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoCompositeShape.h>

#include "AliConst.h"
#include "AliMagF.h"
#include "AliPIPEvGEO.h"
#include "AliRun.h"
#include "AliLog.h"
 
ClassImp(AliPIPEvGEO)
 
//_____________________________________________________________________________
AliPIPEvGEO::AliPIPEvGEO()
{
// Constructor
}

//_____________________________________________________________________________
AliPIPEvGEO::AliPIPEvGEO(const char *name, const char *title)
  : AliPIPE(name,title)
{
// Constructor
}

 
//___________________________________________
void AliPIPEvGEO::CreateGeometry()
{
//Begin_Html
/*
<img src="picts/pipe.gif">
*/
//End_Html


//Begin_Html
/*
<img src="picts/tree_pipe.gif">
*/
//End_Html
//
//
//  The ALICE central beam-pipe according to drawing         LHCVC2C_0001 
//  Drawings of sub-elements:
//  
//  Pos 7 - Minimised Flange:                                LHCVFX_P0025
//  Pos 6 - Standard Flange:                                 STDVFUHV0009
//  Pos 8 - Bellow:                                          LHCVBX__0001
//
//  Absolute z-coordinates -82.0 - 400.0 cm 
//  Total length:                                          482.0 cm
//  It consists of 3 main parts:
//  CP/1 The central Be pipe:                              405.0 cm 
//  CP/2 The flange on the non-absorber side:               36.5 cm  
//  CP/3 The double-bellow and flange on the absorber side: 40.5 cm 
//
//
    AliDebug(1,"Create PIPEv0 geometry");
    Int_t *idtmed = fIdtmed->GetArray();
    Float_t dz;
    Float_t cpcon[39];
    Float_t cpTube[3];
    Int_t   idrotm[2099];  
//
// Rotation Matrices
//
// Rotation by 180 deg
    AliMatrix(idrotm[2012],90.,180., 90., 90.,180.,  0.);

//
    const Float_t kCPz0      = -400.0;
    const Float_t kCP1Length =  405.0;    
    const Float_t kCP2Length =   36.6;
    const Float_t kCP3Length =   40.5;
    
    const Float_t kCP2pos    = kCPz0 + kCP2Length / 2.;
    const Float_t kCP1pos    = kCPz0 + kCP2Length + kCP1Length / 2.;
    const Float_t kCP3pos    = kCPz0 + kCP2Length + kCP1Length + kCP3Length/2.;
///////////////////
//      CP/1     //
///////////////////
//  Inner and outer radii of the Be-section [Pos 1]
    const Float_t kCP1BeRi = 2.90;
    const Float_t kCP1BeRo = 2.98;
//
// Be-Stainless Steel adaptor tube [Pos 2] at both ends of the Be-section. Length 5 cm
    const Float_t kCP1BeStAdaptorLength = 5.00;
//
// Bulge of the Be-Stainless Steel adaptor Tube [Pos 2]
    const Float_t kCP1BeStRo = 3.05;
//
//  Length of bulge [Pos 2]
    const Float_t kCP1BulgeLength = 0.50;
//
//  Distance between bulges [Pos 2]
    const Float_t kCP1BulgeBulgeDistance = 1.00;
//
// Length of Be-pipe
    const Float_t kCP1BeLength =  kCP1Length - 2. *  kCP1BeStAdaptorLength;
//
// CP/1 Mother volume 
    cpTube[0] = 0.;
    cpTube[1] = kCP1BeStRo;
    cpTube[2] = kCP1Length / 2.;
    gMC->Gsvolu("Q1MO","TUBE", idtmed[kAir], cpTube, 3);    
    gMC->Gspos("Q1MO", 1, "ALIC", 0., 0., kCP1pos, 0, "ONLY");

//
// CP/1 Be-Section
//    
    cpTube[0] = 0.;
    cpTube[1] = kCP1BeRo;
    cpTube[2] = kCP1BeLength / 2.;
    gMC->Gsvolu("Q1BE","TUBE", idtmed[kBe], cpTube, 3);    
    cpTube[0] = 0.;
    cpTube[1] = kCP1BeRi;
    gMC->Gsvolu("Q1BV","TUBE", idtmed[kVac], cpTube, 3);    
    gMC->Gspos("Q1BV", 1, "Q1BE", 0., 0., 0., 0, "ONLY");
    gMC->Gspos("Q1BE", 1, "Q1MO", 0., 0., 0., 0, "ONLY");
//
//   CP/1 Be-Stainless Steel adaptor tube
//
    cpcon[0 ]  =   0;
    cpcon[1 ]  = 360;
    cpcon[2 ]  =  8;
//  1 First Bulge 
    cpcon[3 ]  = - kCP1BeStAdaptorLength / 2.;
    cpcon[4 ]  = 0.;
    cpcon[5 ]  = kCP1BeStRo;
//  2
    cpcon[6 ]  = cpcon[3] + kCP1BulgeLength;
    cpcon[7 ]  = 0.;
    cpcon[8 ]  = kCP1BeStRo;
//  3 
    cpcon[9 ]  = cpcon[6];
    cpcon[10]  = 0.;
    cpcon[11]  = kCP1BeRo;
//  4 Between the bulges
    cpcon[12]  = cpcon[9] + kCP1BulgeBulgeDistance;
    cpcon[13]  = 0.;
    cpcon[14]  = kCP1BeRo;
//  5 
    cpcon[15]  = cpcon[12];
    cpcon[16]  = 0.;
    cpcon[17]  = kCP1BeStRo;
//  6 Second bulge 
    cpcon[18]  = cpcon[15] + kCP1BulgeLength;
    cpcon[19]  = 0.;
    cpcon[20]  = kCP1BeStRo;
//  7 
    cpcon[21]  = cpcon[18] + kCP1BulgeLength;
    cpcon[22]  = 0.;
    cpcon[23]  = kCP1BeRo;
// 8 Straight piece
    cpcon[24]  = kCP1BeStAdaptorLength / 2.;
    cpcon[25]  = 0.;
    cpcon[26]  = kCP1BeRo;
//
    gMC->Gsvolu("Q1AT","PCON", idtmed[kInox], cpcon, 27);    
// 
//  Vacuum
    cpTube[0] =  0.;
    cpTube[1] = kCP1BeRi;
    cpTube[2] = -1.;
    gMC->Gsvolu("Q1AV","TUBE", idtmed[kVac], cpTube, 3);
    gMC->Gspos("Q1AV", 1, "Q1AT", 0., 0., 0., 0, "ONLY");
//  Position adaptor tube at both ends
    dz = kCP1Length / 2. -  kCP1BeStAdaptorLength / 2.;
    gMC->Gspos("Q1AT", 1, "Q1MO", 0., 0.,  -dz, 0,            "ONLY");
    gMC->Gspos("Q1AT", 2, "Q1MO", 0., 0.,   dz, idrotm[2012], "ONLY");
//
///////////////////
//      CP/2     //
///////////////////
//
// Fixed Point tube [Pos 5]
//
// Inner and outer radii of the Stainless Steel pipe    
    const Float_t kCP2StRi = 2.90;
    const Float_t kCP2StRo = 2.98;
//  
// Transition to central Be-pipe (Bulge)   
// Length
    const Float_t kCP2BulgeLength = 0.80;
//     
// Bulge outer radius
    const Float_t kCP2BulgeRo = 3.05;
//
// Fixed Point at z = 391.7 (IP)
//
// Position of fixed point
    const Float_t kCP2FixedPointZ = 8.30;
//
// Outer radius of fixed point
    const Float_t kCP2FixedPointRo = 3.50;
//
// Length of fixed point
    const Float_t kCP2FixedPointLength = 0.60;
//
// Fixed Flange [Pos 6]    
//
// Fixed flange outer radius
    const Float_t kCP2FixedFlangeRo = 7.60;
//
// Fixed flange inner radius
    const Float_t kCP2FixedFlangeRi = 3.00;
// Fixed flange inner radius bulge
    const Float_t kCP2FixedFlangeBulgeRi = 2.90;
// Fixed flange lengths of sections at inner radius
    const Float_t  kCP2FixedFlangeRecessLengths[3] ={1., 0.08, 0.9};
// Fixed flange length
    const Float_t kCP2FixedFlangeLength = 1.98;
//
// Fixed flange bulge
// Outer radius
     const Float_t kCP2FixedFlangeBulgeRo = 3.00;
//
// Length    
     const Float_t kCP2FixedFlangeBulgeLength = 2.00;
//
// CP/2 Mother Volume
//
    cpcon[0 ]  =   0;
    cpcon[1 ]  = 360;
    cpcon[2 ]  =  8;
//  1 Flange 
    cpcon[3 ]  = - kCP2Length / 2.;
    cpcon[4 ]  = 0.;
    cpcon[5 ]  = kCP2FixedFlangeRo;
//  2
    cpcon[6 ]  = cpcon[3] + kCP2FixedFlangeLength;
    cpcon[7 ]  = 0.;
    cpcon[8 ]  = kCP2FixedFlangeRo;
//  3 Straight section between Flange and Fixed Point
    cpcon[9 ]  = cpcon[6];
    cpcon[10]  = 0.;
    cpcon[11]  = kCP2FixedFlangeBulgeRo;
//  4 
    cpcon[12]  = cpcon[3] + kCP2FixedPointZ - kCP2FixedPointLength / 2.;
    cpcon[13]  = 0.;
    cpcon[14]  = kCP2FixedFlangeBulgeRo;
//  5 Fixed Point
    cpcon[15]  = cpcon[12];
    cpcon[16]  = 0.;
    cpcon[17]  = kCP2FixedPointRo;
//  6 
    cpcon[18]  = cpcon[15] + kCP2FixedPointLength;
    cpcon[19]  = 0.;
    cpcon[20]  = kCP2FixedPointRo;
//  7 Straight section between Fixed Point and transition bulge
    cpcon[21]  = cpcon[18];
    cpcon[22]  = 0.;
    cpcon[23]  = kCP2BulgeRo;
//  8
    cpcon[24]  = kCP2Length / 2.;
    cpcon[25]  = 0.;
    cpcon[26]  = kCP2BulgeRo;
    gMC->Gsvolu("Q2MO","PCON", idtmed[kAir], cpcon, 27);
    dz =  kCP2pos;
    gMC->Gspos("Q2MO", 1, "ALIC", 0., 0., dz, 0, "ONLY");
//
//  CP/2 Fixed Flange [Pos 6] 
//
    cpcon[0 ]  =   0;
    cpcon[1 ]  = 360;
    cpcon[2 ]  =   4;
//  1 
    cpcon[3 ]  = - kCP2FixedFlangeLength / 2.;
    cpcon[4 ]  = 0.;
    cpcon[5 ]  = kCP2FixedFlangeRo;
//  2 
    cpcon[6 ]  = cpcon[3] + kCP2FixedFlangeRecessLengths[0] + kCP2FixedFlangeRecessLengths[1];
    cpcon[7 ]  = 0.;
    cpcon[8 ]  = kCP2FixedFlangeRo;
//  3 
    cpcon[9 ]  = cpcon[6];
    cpcon[10]  = kCP2FixedFlangeRi;
    cpcon[11]  = kCP2FixedFlangeRo;
//  4 
    cpcon[12]  = - cpcon[3];
    cpcon[13]  = kCP2FixedFlangeRi;
    cpcon[14]  = kCP2FixedFlangeRo;
    gMC->Gsvolu("Q2FL","PCON", idtmed[kInox], cpcon, 15);    
//
//  Vacuum to be placed into Fixed Flange
    cpcon[0 ]  =   0;
    cpcon[1 ]  = 360;
    cpcon[2 ]  =   4;
//  1 
    cpcon[3 ]  = - (kCP2FixedFlangeRecessLengths[0] + kCP2FixedFlangeRecessLengths[1]) / 2.;
    cpcon[4 ]  = 0.;
    cpcon[5 ]  = kCP2FixedFlangeRi;
//  2 
    cpcon[6 ]  = cpcon[3] + kCP2FixedFlangeRecessLengths[0];
    cpcon[7 ]  = 0.;
    cpcon[8 ]  = kCP2FixedFlangeRi; 
//  3 Bulge
    cpcon[9 ]  = cpcon[6];
    cpcon[10]  = 0.;
    cpcon[11]  = kCP2FixedFlangeBulgeRi;
//  4 
    cpcon[12]  = -cpcon[3];
    cpcon[13]  = 0.;
    cpcon[14]  = kCP2FixedFlangeBulgeRi;
    gMC->Gsvolu("Q2V1","PCON", idtmed[kVac], cpcon, 15);   
    dz =  - kCP2FixedFlangeLength / 2. - cpcon[3];
    gMC->Gspos("Q2V1", 1, "Q2FL", 0., 0., dz, 0, "ONLY");
// 
    dz =  - kCP2Length / 2. +  kCP2FixedFlangeLength / 2.;
    gMC->Gspos("Q2FL", 1, "Q2MO", 0., 0., dz, 0, "ONLY");
//
//  CP/2 Beam pipe with fixed point and transition bulges
//    
    cpcon[0 ]  =   0;
    cpcon[1 ]  = 360;
    cpcon[2 ]  =  10;
//  1 Bulge at transition to flange 
    cpcon[3 ]  = - (kCP2Length -  kCP2FixedFlangeRecessLengths[0] - kCP2FixedFlangeRecessLengths[1]) / 2.;
    cpcon[4 ]  = 0.;
    cpcon[5 ]  = kCP2FixedFlangeBulgeRo;
//  2
    cpcon[6 ]  = cpcon[3] +  kCP2FixedFlangeBulgeLength;
    cpcon[7 ]  = 0.;
    cpcon[8 ]  = kCP2FixedFlangeBulgeRo;
//  3 Straight section between Bulge and Fixed Point
    cpcon[9 ]  = cpcon[6];
    cpcon[10]  = 0.;
    cpcon[11]  = kCP2StRo;
//  4 
    cpcon[12]  = cpcon[3] + (kCP2FixedPointZ - kCP2FixedFlangeRecessLengths[0] - kCP2FixedFlangeRecessLengths[1])  - kCP2FixedPointLength / 2.;
    cpcon[13]  = 0.;
    cpcon[14]  = kCP2StRo;
//  5 Fixed Point
    cpcon[15]  = cpcon[12];
    cpcon[16]  = 0.;
    cpcon[17]  = kCP2FixedPointRo;
//  6 
    cpcon[18]  = cpcon[15] + kCP2FixedPointLength;
    cpcon[19]  = 0.;
    cpcon[20]  = kCP2FixedPointRo;
//  7 Straight section between Fixed Point and transition bulge
    cpcon[21]  = cpcon[18];
    cpcon[22]  = 0.;
    cpcon[23]  = kCP2StRo;
//  8
    cpcon[24]  = - cpcon[3] - kCP2BulgeLength;
    cpcon[25]  = 0.;
    cpcon[26]  = kCP2StRo;
// 9 Bulge at transition to Be pipe
    cpcon[27]  = cpcon[24];
    cpcon[28]  = 0.;
    cpcon[29]  = kCP2BulgeRo;
// 10 
    cpcon[30]  = - cpcon[3];
    cpcon[31]  = 0.;
    cpcon[32]  = kCP2BulgeRo;
    gMC->Gsvolu("Q2PI","PCON", idtmed[kInox], cpcon, 33);    
//
//  Vacuum to be place into CP/2 beam pipe
    cpTube[0] = 0.;
    cpTube[1] = kCP2StRi;
    cpTube[2] = -1.;
    gMC->Gsvolu("Q2V2","TUBE", idtmed[kVac], cpTube, 3);   

    gMC->Gspos("Q2V2", 1, "Q2PI", 0., 0., 0., 0, "ONLY");    
    dz = (kCP2FixedFlangeRecessLengths[0] + kCP2FixedFlangeRecessLengths[1]) / 2.;
    gMC->Gspos("Q2PI", 1, "Q2MO", 0., 0., dz, 0, "ONLY");    
//
///////////////////
//      CP/3     //
///////////////////
//
// Adaptor tube [Pos 4]
// 
// Adaptor tube length 
    const Float_t  kCP3AdaptorTubeLength = 5.50;
//
// Inner and outer radii
     const Float_t kCP3AdaptorTubeRi = 2.92;
     const Float_t kCP3AdaptorTubeRo = 3.00;
//
// Bulge at transition point
// Inner and outer radii
     const Float_t kCP3AdaptorTubeBulgeRi = 2.90;
     const Float_t kCP3AdaptorTubeBulgeRo = 3.05;    
//
// Length of bulge
    const Float_t  kCP3AdaptorTubeBulgeLength = 0.80;
//
// Bellow [Pos 8]
//
//  Total length    
    const Float_t kCP3BellowLength = 13.00;
//  Outer Radius
    const Float_t kCP3BellowRo = 3.6;
//  Inner Radius 
    const Float_t kCP3BellowRi = 2.8;
//  Number of plies
    const Int_t   kCP3NumberOfPlies = 18;
//  Length of undulated region
    const Float_t kCP3BellowUndulatedLength  = 8.30; 
//  Plie thickness
    const Float_t kCP3PlieThickness = 0.02;   
//  Connection Plie radies (at transition been undulated region and beam pipe)
    const Float_t kCP3ConnectionPlieR = 0.21;
//  Plie radius
//  const Float_t kCP3PlieR = 0.118286;
    const Float_t kCP3PlieR = 
	(kCP3BellowUndulatedLength - 4. *  kCP3ConnectionPlieR + 2. * kCP3PlieThickness + (2. *  kCP3NumberOfPlies - 2.) * kCP3PlieThickness) 
	/ (4. * kCP3NumberOfPlies - 2.);
//  Length of connection pipe
    const Float_t kCP3BellowConnectionLength = 2.35;
//
//  Tube between bellows [Pos 3]  
//    
//  Length of tube
    const Float_t kCP3TubeLength = 4.00;
//
//  Minimised fixed flange [Pos 7]
//  
//  Length of flange connection tube
    const Float_t kCP3FlangeConnectorLength = 5.0 - 1.4;
//  Length of Flange
    const Float_t kCP3FlangeLength = 1.40;
//  Outer radius    
    const Float_t kCP3FlangeRo     = 4.30;

//
// CP/3 Mother volume
//
    cpcon[0 ]  =   0;
    cpcon[1 ]  = 360;
    cpcon[2 ]  =  12;
//  1 From transition to first bellow
    cpcon[3 ]  = - kCP3Length / 2.;
    cpcon[4 ]  = 0.;
    cpcon[5 ]  = kCP3AdaptorTubeBulgeRo;
//  2
    cpcon[6 ]  = cpcon[3] + kCP3BellowConnectionLength + kCP3AdaptorTubeLength;
    cpcon[7 ]  = 0.;
    cpcon[8 ]  = kCP3AdaptorTubeBulgeRo;
//  3 First Bellow
    cpcon[9 ]  = cpcon[6];
    cpcon[10]  = 0.;
    cpcon[11]  = kCP3BellowRo;
//  4 
    cpcon[12]  = cpcon[9] +  kCP3BellowUndulatedLength;
    cpcon[13]  = 0.;
    cpcon[14]  = kCP3BellowRo;
//  5 Connection between the two bellows
    cpcon[15]  = cpcon[12];
    cpcon[16]  = 0.;
    cpcon[17]  = kCP3AdaptorTubeBulgeRo;
//  6 
    cpcon[18]  = cpcon[15] + 2. * kCP3BellowConnectionLength + kCP3TubeLength;
    cpcon[19]  = 0.;
    cpcon[20]  = kCP3AdaptorTubeBulgeRo;
//  7 Second bellow
    cpcon[21]  = cpcon[18];
    cpcon[22]  = 0.;
    cpcon[23]  = kCP3BellowRo;
//  8
    cpcon[24]  = cpcon[21] + kCP3BellowUndulatedLength;
    cpcon[25]  = 0.;
    cpcon[26]  = kCP3BellowRo;
//  9 Pipe between second Bellow and Flange 
    cpcon[27]  = cpcon[24];
    cpcon[28]  = 0.;
    cpcon[29]  = kCP3AdaptorTubeBulgeRo;
//  10 
    cpcon[30]  = cpcon[27] + kCP3BellowConnectionLength +  kCP3FlangeConnectorLength;
    cpcon[31]  = 0.;
    cpcon[32]  = kCP3AdaptorTubeBulgeRo;
//  11 Flange 
    cpcon[33]  = cpcon[30];
    cpcon[34]  = 0.;
    cpcon[35]  = kCP3FlangeRo;
//  12 
    cpcon[36]  = - cpcon[3];
    cpcon[37]  = 0.;
    cpcon[38]  = kCP3FlangeRo;
//
    gMC->Gsvolu("Q3MO","PCON", idtmed[kAir], cpcon, 39);
    dz = kCP3pos;
    gMC->Gspos("Q3MO", 1, "ALIC", 0., 0., dz, 0, "ONLY");
//
// CP/3 Adaptor tube
//
    cpcon[0 ]  =   0;
    cpcon[1 ]  = 360;
    cpcon[2 ]  =  4;
//  1 Bulge at transition
    cpcon[3 ]  = - kCP3AdaptorTubeLength / 2.;
    cpcon[4 ]  = 0.;
    cpcon[5 ]  = kCP3AdaptorTubeBulgeRo;
//  2
    cpcon[6 ]  = cpcon[3] + kCP3AdaptorTubeBulgeLength;
    cpcon[7 ]  = 0.;
    cpcon[8 ]  = kCP3AdaptorTubeBulgeRo;
//  3 Tube
    cpcon[9 ]  = cpcon[6];
    cpcon[10]  = 0.;
    cpcon[11]  = kCP3AdaptorTubeRo;
//  4 
    cpcon[12]  = - cpcon[3];
    cpcon[13]  = 0.;
    cpcon[14]  = kCP3AdaptorTubeRo;
    gMC->Gsvolu("Q3ATO","PCON", idtmed[kVac], cpcon, 15);

//  1 Bulge at transition
    cpcon[4 ]  = kCP3AdaptorTubeBulgeRi;
//  2
    cpcon[7 ]  = kCP3AdaptorTubeBulgeRi;;
//  3 Tube
    cpcon[10]  = kCP3AdaptorTubeRi;
//  4 
    cpcon[13]  = kCP3AdaptorTubeRi;

    gMC->Gsvolu("Q3ATI","PCON", idtmed[kInox], cpcon, 15);
    gMC->Gspos("Q3ATI", 1, "Q3ATO", 0., 0., 0., 0, "ONLY");
    dz = - kCP3Length / 2. +  kCP3AdaptorTubeLength / 2.;
    gMC->Gspos("Q3ATO", 1, "Q3MO", 0., 0., dz, 0, "ONLY");
//
// CP/3 Bellow section
//
//    TGeoMedium* mAir    =  gGeoManager->GetMedium("AIR");
    TGeoMedium* mVacuum =  gGeoManager->GetMedium("VACUUM");    
    TGeoMedium* mSteel  =  gGeoManager->GetMedium("INOX");        
//
//  Upper part of the undulation
    TGeoTorus* plieTorusUO =  new TGeoTorus(kCP3BellowRo - kCP3PlieR, 0. , kCP3PlieR);
    plieTorusUO->SetName("TorusUO");
    TGeoTorus* plieTorusUI =  new TGeoTorus(kCP3BellowRo - kCP3PlieR, kCP3PlieR - kCP3PlieThickness, kCP3PlieR);
    plieTorusUI->SetName("TorusUI");
    TGeoTube*  plieTubeU   =  new TGeoTube (kCP3BellowRo - kCP3PlieR, kCP3BellowRo, kCP3PlieR / 2.);
    plieTubeU->SetName("TubeU");
    
    TGeoCompositeShape*  upperPlieO = new TGeoCompositeShape("upperPlieO", "TorusUO*TubeU");
    TGeoCompositeShape*  upperPlieI = new TGeoCompositeShape("upperPlieI", "TorusUI*TubeU");
 
    TGeoVolume* wiggleUO = new TGeoVolume("Q3WUO", upperPlieO, mVacuum);
    TGeoVolume* wiggleUI = new TGeoVolume("Q3WUI", upperPlieI, mSteel);
    wiggleUO->AddNode(wiggleUI, 1,  new TGeoTranslation(0., 0., 0.));    
//
// Lower part of the undulation
    TGeoTorus* plieTorusLO =  new TGeoTorus(kCP3BellowRi + kCP3PlieR, 0. , kCP3PlieR);
    plieTorusLO->SetName("TorusLO");
    TGeoTorus* plieTorusLI =  new TGeoTorus(kCP3BellowRi + kCP3PlieR, kCP3PlieR - kCP3PlieThickness, kCP3PlieR);
    plieTorusLI->SetName("TorusLI");
    TGeoTube*  plieTubeL   =  new TGeoTube (kCP3BellowRi, kCP3BellowRi + kCP3PlieR, kCP3PlieR);
    plieTubeL->SetName("TubeL");

    TGeoCompositeShape*  lowerPlieO = new TGeoCompositeShape("lowerPlieO", "TorusLO*TubeL");
    TGeoCompositeShape*  lowerPlieI = new TGeoCompositeShape("lowerPlieI", "TorusLI*TubeL");

    TGeoVolume* wiggleLO = new TGeoVolume("Q3WLO", lowerPlieO, mVacuum);
    TGeoVolume* wiggleLI = new TGeoVolume("Q3WLI", lowerPlieI, mSteel);
    wiggleLO->AddNode(wiggleLI, 1,  new TGeoTranslation(0., 0., 0.));    

//
// Connection between upper and lower part of undulation
    TGeoVolume* wiggleC1 = new TGeoVolume("Q3WCO1",  
					  new TGeoTube(kCP3BellowRi + kCP3PlieR, kCP3BellowRo - kCP3PlieR, kCP3PlieThickness / 2.),
					  mSteel);
    TGeoVolume* wiggleC2 = new TGeoVolume("Q3WCO2",  
					  new TGeoTube(kCP3BellowRi + kCP3ConnectionPlieR, kCP3BellowRo - kCP3PlieR, kCP3PlieThickness / 2.),
					  mSteel);
//
// Conncetion between undulated section and beam pipe
    TGeoTorus* plieTorusCO =  new TGeoTorus(kCP3BellowRi + kCP3ConnectionPlieR, 0. , kCP3ConnectionPlieR);
    plieTorusCO->SetName("TorusCO");
    TGeoTorus* plieTorusCI =  new TGeoTorus(kCP3BellowRi + kCP3ConnectionPlieR, kCP3ConnectionPlieR - kCP3PlieThickness, kCP3ConnectionPlieR);
    plieTorusCI->SetName("TorusCI");
    TGeoTube*  plieTubeC   =  new TGeoTube (kCP3BellowRi, kCP3BellowRi + kCP3ConnectionPlieR, kCP3ConnectionPlieR);
    plieTubeC->SetName("TubeC");

    TGeoCompositeShape*  connectionPlieO = new TGeoCompositeShape("connectionPlieO", "TorusCO*TubeC");
    TGeoCompositeShape*  connectionPlieI = new TGeoCompositeShape("connectionPlieI", "TorusCI*TubeC");

    TGeoVolume* connectionPO = new TGeoVolume("Q3CPO", connectionPlieO, mVacuum);
    TGeoVolume* connectionPI = new TGeoVolume("Q3CPI", connectionPlieI, mSteel);
    connectionPO->AddNode(connectionPI, 1,  new TGeoTranslation(0., 0., 0.));    
//
// Connecting pipes
    TGeoVolume* connectionPipeO = new TGeoVolume("Q3BECO",  
						 new TGeoTube(0., kCP3AdaptorTubeRo, kCP3BellowConnectionLength / 2.),
						 mVacuum);
    TGeoVolume* connectionPipeI = new TGeoVolume("Q3BECI",  
						 new TGeoTube(kCP3AdaptorTubeRi, kCP3AdaptorTubeRo, kCP3BellowConnectionLength / 2.),
						 mSteel);
  
    connectionPipeO->AddNode(connectionPipeI, 1,  new TGeoTranslation(0., 0., 0.));
    
//
// Bellow mother
    TGeoPcon* bellowMotherPC = new TGeoPcon(0., 360., 6);
    dz =  - kCP3BellowLength / 2;
    bellowMotherPC->DefineSection(0, dz, 0.,  kCP3AdaptorTubeRo);
    dz +=  kCP3BellowConnectionLength;
    bellowMotherPC->DefineSection(1, dz, 0.,  kCP3AdaptorTubeRo);
    bellowMotherPC->DefineSection(2, dz, 0.,  kCP3BellowRo);
    dz =  kCP3BellowLength /2. -  kCP3BellowConnectionLength;;
    bellowMotherPC->DefineSection(3, dz, 0.,  kCP3BellowRo);
    bellowMotherPC->DefineSection(4, dz, 0.,  kCP3AdaptorTubeRo);
    dz +=  kCP3BellowConnectionLength;
    bellowMotherPC->DefineSection(5, dz, 0.,  kCP3AdaptorTubeRo);

    TGeoVolume* bellowMother = new TGeoVolume("Q3BeMO", bellowMotherPC, mVacuum);
//
// Add undulations
    Float_t z0 =  - kCP3BellowLength / 2. +  kCP3BellowConnectionLength + 2. * kCP3ConnectionPlieR - kCP3PlieThickness;
    Float_t zsh  = 4. *  kCP3PlieR -  2. * kCP3PlieThickness;
    for (Int_t iw = 0; iw < 18; iw++) {
	Float_t zpos =  z0 + iw * zsh;	
	if (iw > 0) 
	    bellowMother->AddNode(wiggleC1,  iw + 1 , new TGeoTranslation(0., 0., zpos + kCP3PlieThickness / 2.));	
	else
	    bellowMother->AddNode(wiggleC2,  iw + 1 , new TGeoTranslation(0., 0., zpos + kCP3PlieThickness / 2.));	

	zpos += kCP3PlieR;
	bellowMother->AddNode(wiggleUO, iw + 1,  new TGeoTranslation(0., 0., zpos));	

	zpos += kCP3PlieR;
	if (iw < 17) 
	    bellowMother->AddNode(wiggleC1,  iw + 19, new TGeoTranslation(0., 0., zpos - kCP3PlieThickness / 2.));
	else
	    bellowMother->AddNode(wiggleC2,  iw + 19, new TGeoTranslation(0., 0., zpos - kCP3PlieThickness / 2.));

	if (iw < 17) {
	    zpos += kCP3PlieR;
	    bellowMother->AddNode(wiggleLO, iw + 1, new TGeoTranslation(0., 0., zpos -  kCP3PlieThickness));
	}
    }
//
// Add connecting undulation between bellow and connecting pipe
    dz = - kCP3BellowUndulatedLength / 2. + kCP3ConnectionPlieR;
    bellowMother->AddNode(connectionPO, 1,  new TGeoTranslation(0., 0.,  dz));
    bellowMother->AddNode(connectionPO, 2,  new TGeoTranslation(0., 0., -dz));
//
// Add connecting pipe
    dz =  - kCP3BellowLength / 2. +  kCP3BellowConnectionLength / 2.;
    bellowMother->AddNode(connectionPipeO, 1,  new TGeoTranslation(0., 0.,   dz));
    bellowMother->AddNode(connectionPipeO, 2,  new TGeoTranslation(0., 0.,  -dz));
//
// Add bellow to CP/3 mother    
    TGeoVolume* mother = gGeoManager->GetVolume("Q3MO");
    dz = - kCP3Length / 2. +  kCP3AdaptorTubeLength +  kCP3BellowLength / 2.;
    mother->AddNode(bellowMother, 1,  new TGeoTranslation(0., 0., dz));
    dz += (kCP3BellowLength +  kCP3TubeLength);
    mother->AddNode(bellowMother, 2,  new TGeoTranslation(0., 0., dz));
//
// Beam pipe section between bellows
//
    cpTube[0] = 0.;
    cpTube[1] = kCP3AdaptorTubeRo;
    cpTube[2] = kCP3TubeLength / 2.;
    gMC->Gsvolu("Q3BCO","TUBE", idtmed[kVac], cpTube, 3);

    cpTube[0] = kCP3AdaptorTubeRi;
    cpTube[1] = kCP3AdaptorTubeRo;
    cpTube[2] = kCP3TubeLength / 2.;
    gMC->Gsvolu("Q3BCI","TUBE", idtmed[kInox], cpTube, 3);
  
    gMC->Gspos("Q3BCI", 1, "Q3BCO", 0., 0., 0., 0, "ONLY");
    dz = - kCP3Length / 2. +   kCP3AdaptorTubeLength +  kCP3BellowLength +  kCP3TubeLength / 2.;
    gMC->Gspos("Q3BCO", 1, "Q3MO", 0., 0., dz, 0, "ONLY");
    

// CP3 Minimised Flange
//
    cpcon[0 ]  =   0;
    cpcon[1 ]  = 360;
    cpcon[2 ]  =  4;
//  1 Connection Tube
    cpcon[3 ]  = - (kCP3FlangeConnectorLength + kCP3FlangeLength) / 2.;
    cpcon[4 ]  = 0.;
    cpcon[5 ]  = kCP3AdaptorTubeRo;
//  2
    cpcon[6 ]  = cpcon[3] + kCP3FlangeConnectorLength;
    cpcon[7 ]  = 0.;
    cpcon[8 ]  = kCP3AdaptorTubeRo;
//  3 Flange
    cpcon[9 ]  = cpcon[6];
    cpcon[10]  = 0.;
    cpcon[11]  = kCP3FlangeRo;
//  4 
    cpcon[12]  = - cpcon[3];
    cpcon[13]  = 0.;
    cpcon[14]  = kCP3FlangeRo;
    gMC->Gsvolu("Q3MFO","PCON", idtmed[kVac], cpcon, 15);

    cpcon[4 ]  = cpcon[7 ] = cpcon [10] = cpcon[13] = kCP3AdaptorTubeRi;
    gMC->Gsvolu("Q3MFI","PCON", idtmed[kInox], cpcon, 15);

    gMC->Gspos("Q3MFI", 1, "Q3MFO", 0., 0., 0., 0, "ONLY");
    dz =  kCP3Length / 2. - (kCP3FlangeConnectorLength + kCP3FlangeLength) / 2.;
    gMC->Gspos("Q3MFO", 1, "Q3MO", 0., 0., dz, 0, "ONLY");
}



//___________________________________________
void AliPIPEvGEO::CreateMaterials()
{
  //
  // Define materials for beam pipe
  //

  AliDebugClass(1,"Create PIPEvGEO materials");
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  // Steel (Inox)  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  // AlBe - alloy 
  Float_t aAlBe[2] = { 26.98, 9.01};
  Float_t zAlBe[2] = { 13.00, 4.00};
  Float_t wAlBe[2] = { 0.4, 0.6};
  //
  // Polyamid
  Float_t aPA[4] = {16., 14., 12.,  1.};
  Float_t zPA[4] = { 8.,  7.,  6.,  1.};
  Float_t wPA[4] = { 1.,  1.,  6., 11.};
  //
  // Air 
  //
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  Float_t dAir1 = 1.20479E-10;

  //
  //     Berillium 
  AliMaterial(5, "BERILLIUM$", 9.01, 4., 1.848, 35.3, 36.7);
  //
  //     Carbon 
  AliMaterial(6,  "CARBON$   ", 12.01, 6., 2.265, 18.8, 49.9);
  //
  //     Aluminum 
  AliMaterial(9,  "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  //
  //     Air 
  AliMixture(15, "AIR$      ", aAir, zAir, dAir, 4, wAir);
  //
  //     Vacuum 
  AliMixture(16, "VACUUM$ ", aAir, zAir, dAir1, 4, wAir);
  //
  //     stainless Steel 
  AliMixture(19, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  //
  //     reduced density steel to approximate pump getter material
  AliMixture(20, "GETTER$", asteel, zsteel, 1.00, 4, wsteel);
  //     Al-Be alloy
  //     
  AliMixture(21, "AlBe$", aAlBe, zAlBe, 2.07, 2, wAlBe);
  //     Polyamid
  //   
  AliMixture(22, "PA$", aPA, zPA, 1.14, -4, wPA);
  //

  // **************** 
  //     Defines tracking media parameters. 
  //
  Float_t epsil  = .001;    // Tracking precision, 
  Float_t stemax = -0.01;   // Maximum displacement for multiple scat 
  Float_t tmaxfd = -20.;    // Maximum angle due to field deflection 
  Float_t deemax = -.3;     // Maximum fractional energy loss, DLS 
  Float_t stmin  = -.8;
  // *************** 
  //
  //    Beryllium 
  
  AliMedium(5, "BE",       5, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

  //    Carbon 
  AliMedium(6, "C",        6, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Aluminum 
  AliMedium(9, "ALU",      9, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Air 
  AliMedium(15, "AIR",    15, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Vacuum 
  AliMedium(16, "VACUUM", 16, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Steel 
  AliMedium(19, "INOX",   19, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Getter 
  AliMedium(20, "GETTER", 20, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //   AlBe - Aloy 
  AliMedium(21, "AlBe"  , 21, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //   Polyamid
  AliMedium(22, "PA"  ,   22, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

}










