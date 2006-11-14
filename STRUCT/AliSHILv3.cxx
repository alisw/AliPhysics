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

/* $Id $ */

//-------------------------------------------------------------------------
// The small angle absorber SAA (beam shield)
// Author: A.Morsch
// andreas.morsch@cern.ch
//-------------------------------------------------------------------------

#include <TVirtualMC.h>
#include <TArrayI.h>
#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>
#include <TGeoBBox.h>
#include <TGeoPgon.h>
#include <TGeoCone.h>

#include "AliSHILv3.h"
#include "AliConst.h"
#include "AliLog.h"

ClassImp(AliSHILv3)
 
//_____________________________________________________________________________
AliSHILv3::AliSHILv3()
{
  //
  // Default constructor for muon shield
  //
}
 
//_____________________________________________________________________________
AliSHILv3::AliSHILv3(const char *name, const char *title)
  : AliSHIL(name,title)
{
  //
  // Standard constructor for muon shield
  //
}
 
//_____________________________________________________________________________
void AliSHILv3::CreateGeometry()
{
//
// The geometry of the small angle absorber "Beam Shield"
//
    Float_t dz, dr, z, rmax;
//
// The top volume
//
    TGeoVolume* top = gGeoManager->GetVolume("ALIC");
    
//
// Translations
//
    TGeoTranslation* vec0 = new TGeoTranslation(0., 0., 0.);
//  Rotations
    TGeoRotation* rot000 = new TGeoRotation("rot000",  90.,   0., 90.,  90., 0., 0.);
    TGeoRotation* rot090 = new TGeoRotation("rot090",  90.,  90., 90., 180., 0., 0.);
    TGeoRotation* rot180 = new TGeoRotation("rot180",  90., 180., 90., 270., 0., 0.);
    TGeoRotation* rot270 = new TGeoRotation("rot270",  90., 270., 90.,   0., 0., 0.);
//
// Media
//
    TGeoMedium* kMedNiW     = gGeoManager->GetMedium("SHIL_Ni/W0");
    TGeoMedium* kMedNiWsh   = gGeoManager->GetMedium("SHIL_Ni/W3");
    TGeoMedium* kMedSteel   = gGeoManager->GetMedium("SHIL_ST_C0");
    TGeoMedium* kMedSteelSh = gGeoManager->GetMedium("SHIL_ST_C3");
    TGeoMedium* kMedAir     = gGeoManager->GetMedium("SHIL_AIR_C0");
    TGeoMedium* kMedAirMu   = gGeoManager->GetMedium("SHIL_AIR_MUON");
    TGeoMedium* kMedPb      = gGeoManager->GetMedium("SHIL_PB_C0");
    TGeoMedium* kMedPbSh    = gGeoManager->GetMedium("SHIL_PB_C2");
    TGeoMedium* kMedConc    = gGeoManager->GetMedium("SHIL_CC_C0");
//
    const Float_t kDegRad = TMath::Pi() / 180.;
    const Float_t kAngle02   = TMath::Tan( 2.00   * kDegRad);   
    const Float_t kAngle0071 = TMath::Tan( 0.71   * kDegRad);   


///////////////////////////////////
//    FA Tungsten Tail           //
//    Drawing ALIP2A__0049       //
//    Drawing ALIP2A__0111       //
///////////////////////////////////
//
//    The tail as built is shorter than in drawing ALIP2A__0049. 
//    The CDD data base has to be updated !
//
//    Inner radius at the entrance of the flange
      Float_t rInFaWTail1  = 13.98/2.;
//    Outer radius at the entrance of the flange
      Float_t rOuFaWTail1  = 52.00/2.;
//    Outer radius at the end of the section inside the FA
      Float_t rOuFaWTail2  = 35.27/2.;
//    Length of the Flange section inside the FA
      Float_t dzFaWTail1   =  6.00;
//    Length of the Flange section ouside the FA
      Float_t dzFaWTail2   = 12.70;
//    Inner radius at the end of the section inside the FA 
      Float_t rInFaWTail2  = rInFaWTail1 +  dzFaWTail1 * kAngle0071;
//    Inner radius at the end of the flange
      Float_t rInFaWTail3  = rInFaWTail2 +  dzFaWTail2 * kAngle0071;
//    Outer radius at the end of the flange
      Float_t rOuFaWTail3  = rOuFaWTail2 +  dzFaWTail2 * kAngle02;
//    Outer radius of the recess for station 1
      Float_t rOuFaWTailR  = 30.8/2.;
//    Length of the recess
      Float_t dzFaWTailR   = 36.00;
//    Inner radiues at the end of the recess      
      Float_t rInFaWTail4  =  rInFaWTail3 +  dzFaWTailR * kAngle0071;
//    Outer radius at the end of the recess      
      Float_t rOuFaWTail4  =  rOuFaWTail3 +  dzFaWTailR * kAngle02;
//    Inner radius of the straight section
      Float_t rInFaWTailS  = 22.30/2.;
//    Length of the bulge
      Float_t dzFaWTailB   = 13.0;
//    Outer radius at the end of the bulge
      Float_t rOuFaWTailB  =  rOuFaWTail4 +  dzFaWTailB * kAngle02;
//    Outer radius at the end of the tail 
      Float_t rOuFaWTailE  = 31.6/2.;
//    Total length of the tail
      Float_t dzFaWTail    = 70.7;

      TGeoPcon* shFaWTail = new TGeoPcon(0., 360., 10);
      z    = 0.;
//    Flange section inside FA
      shFaWTail->DefineSection(0,  z, rInFaWTail1, rOuFaWTail1);
      z += dzFaWTail1;
      shFaWTail->DefineSection(1,  z, rInFaWTail2, rOuFaWTail1);      
      shFaWTail->DefineSection(2,  z, rInFaWTail2, rOuFaWTail2);      
//    Flange section outside FA
      z +=  dzFaWTail2;
      shFaWTail->DefineSection(3,  z, rInFaWTail3, rOuFaWTail3);      
      shFaWTail->DefineSection(4,  z, rInFaWTail3, rOuFaWTailR);
//    Recess Station 1
      z +=  dzFaWTailR;
      shFaWTail->DefineSection(5,  z, rInFaWTail4, rOuFaWTailR); 
      shFaWTail->DefineSection(6,  z, rInFaWTailS, rOuFaWTail4);
//    Bulge
      z += dzFaWTailB;
      shFaWTail->DefineSection(7,  z, rInFaWTailS, rOuFaWTailB);
      shFaWTail->DefineSection(8,  z, rInFaWTailS, rOuFaWTailE);
//    End
      z =  dzFaWTail;
      shFaWTail->DefineSection(9,  z, rInFaWTailS, rOuFaWTailE);
      
      TGeoVolume* voFaWTail = new TGeoVolume("YFaWTail", shFaWTail, kMedNiW);
//
//    Define an inner region with higher transport cuts
      TGeoPcon* shFaWTailI = new TGeoPcon(0., 360., 4);
      z  = 0.;
      dr = 3.5;
      shFaWTailI->DefineSection(0, z, rInFaWTail1,  rInFaWTail1 + dr);
      z += (dzFaWTail1 + dzFaWTail2 +  dzFaWTailR);
      shFaWTailI->DefineSection(1, z, rInFaWTail4,  rInFaWTail4 + dr);
      shFaWTailI->DefineSection(2, z, rInFaWTailS,  rInFaWTailS + dr);
      z = dzFaWTail;
      shFaWTailI->DefineSection(3, z, rInFaWTailS,  rInFaWTailS + dr);
      TGeoVolume* voFaWTailI = new TGeoVolume("YFaWTailI", shFaWTailI, kMedNiWsh);
      voFaWTail->AddNode(voFaWTailI, 1, vec0);
      
///////////////////////////////////
//                               //
// Recess Station 1              //
// Drawing ALIP2A__0260          //
///////////////////////////////////

///////////////////////////////////
//    FA W-Ring 2                //
//    Drawing ALIP2A__0220       //
///////////////////////////////////
      const Float_t kFaWring2Rinner  = 15.40;
      const Float_t kFaWring2Router  = 18.40;
      const Float_t kFaWring2HWidth  =  3.75;
      const Float_t kFaWring2Cutoffx =  3.35;
      const Float_t kFaWring2Cutoffy =  3.35;
      TGeoTubeSeg* shFaWring2a  = new TGeoTubeSeg(kFaWring2Rinner, kFaWring2Router, kFaWring2HWidth, 0., 90.);
      shFaWring2a->SetName("shFaWring2a");
      TGeoBBox* shFaWring2b  = new TGeoBBox(kFaWring2Router / 2., kFaWring2Router / 2., kFaWring2HWidth);
      shFaWring2b->SetName("shFaWring2b");
      TGeoTranslation* trFaWring2b 
	  = new TGeoTranslation("trFaWring2b", kFaWring2Router / 2. + kFaWring2Cutoffx, kFaWring2Router / 2. + kFaWring2Cutoffy, 0.);
      trFaWring2b->RegisterYourself();
      TGeoCompositeShape*  shFaWring2 = new TGeoCompositeShape("shFaWring2", "(shFaWring2a)*(shFaWring2b:trFaWring2b)");
      TGeoVolume* voFaWring2    = new TGeoVolume("YFA_WRING2", shFaWring2, kMedNiW);

///////////////////////////////////
//    FA W-Ring 3                //
//    Drawing ALIP2A__0219       //
///////////////////////////////////
      const Float_t kFaWring3Rinner  = 15.40;
      const Float_t kFaWring3Router  = 18.40;
      const Float_t kFaWring3HWidth  =  3.75;
      const Float_t kFaWring3Cutoffx =  3.35;
      const Float_t kFaWring3Cutoffy =  3.35;
      TGeoTubeSeg* shFaWring3a  = new TGeoTubeSeg(kFaWring3Rinner, kFaWring3Router, kFaWring3HWidth, 0., 90.);
      shFaWring3a->SetName("shFaWring3a");
      TGeoBBox* shFaWring3b  = new TGeoBBox(kFaWring3Router / 2., kFaWring3Router / 2., kFaWring3HWidth);
      shFaWring3b->SetName("shFaWring3b");
      TGeoTranslation* trFaWring3b 
	  = new TGeoTranslation("trFaWring3b", kFaWring3Router / 2. + kFaWring3Cutoffx, kFaWring3Router / 2. + kFaWring3Cutoffy, 0.);
      trFaWring3b->RegisterYourself();
      TGeoCompositeShape*  shFaWring3 = new TGeoCompositeShape("shFaWring3", "(shFaWring3a)*(shFaWring3b:trFaWring3b)");
      TGeoVolume* voFaWring3    = new TGeoVolume("YFA_WRING3", shFaWring3, kMedNiW);

///////////////////////////////////
//    FA W-Ring 5                //
//    Drawing ALIP2A__0221       //
///////////////////////////////////
      const Float_t kFaWring5Rinner = 15.40;
      const Float_t kFaWring5Router = 18.67;
      const Float_t kFaWring5HWidth =  1.08;
      TGeoVolume* voFaWring5    = new TGeoVolume("YFA_WRING5", new TGeoTube(kFaWring5Rinner, kFaWring5Router, kFaWring5HWidth), kMedNiW);

//
// Position the rings in the assembly 
//      
      TGeoVolumeAssembly* asFaExtraShield = new TGeoVolumeAssembly("YCRE");
// Distance between rings
      const Float_t kFaDWrings = 1.92;
//
      dz = 0.;
      
      dz +=  kFaWring2HWidth;
      asFaExtraShield->AddNode(voFaWring2,    1, new TGeoCombiTrans(0., 0., dz, rot090));
      asFaExtraShield->AddNode(voFaWring2,    2, new TGeoCombiTrans(0., 0., dz, rot270));
      dz +=   kFaWring2HWidth;
      dz +=   kFaDWrings;
      dz +=   kFaWring3HWidth;
      asFaExtraShield->AddNode(voFaWring3,    1, new TGeoCombiTrans(0., 0., dz, rot000));
      asFaExtraShield->AddNode(voFaWring3,    2, new TGeoCombiTrans(0., 0., dz, rot180));
      dz +=   kFaWring3HWidth;   
      dz +=   kFaWring5HWidth;   
      asFaExtraShield->AddNode(voFaWring5,    1, new TGeoTranslation(0., 0., dz));
      dz +=   kFaWring5HWidth;   
      dz +=   kFaWring3HWidth;   
      asFaExtraShield->AddNode(voFaWring3,    3, new TGeoCombiTrans(0., 0., dz, rot090));
      asFaExtraShield->AddNode(voFaWring3,    4, new TGeoCombiTrans(0., 0., dz, rot270));
      dz +=   kFaWring3HWidth;   
      dz +=   kFaDWrings;
      dz +=   kFaWring2HWidth;
      asFaExtraShield->AddNode(voFaWring2,    3, new TGeoCombiTrans(0., 0., dz, rot000));
      asFaExtraShield->AddNode(voFaWring2,    4, new TGeoCombiTrans(0., 0., dz, rot180));
      dz +=   kFaWring2HWidth;

      
///////////////////////////////////////
//                SAA1               //
///////////////////////////////////////


///////////////////////////////////////
//          FA/SAA1  W Joint         //
//          Drawing ALIP2A__0060     //
///////////////////////////////////////

// Length of flange FA side      
      Float_t dzFaSaa1F1  =  2.8;
// Inner radius of flange FA side
      Float_t rInFaSaa1F1 = 32.0/2.;
// Outer radius of flange FA side
      Float_t rOuFaSaa1F1 = 39.5/2.;
// Length of first straight section
      Float_t dzFaSaa1S1  = 18.5 - dzFaSaa1F1;
// Inner radius of first straight section
      Float_t rInFaSaa1S1 = 22.3/2.;
// Length of 45 deg transition region
      Float_t  dzFaSaa1T1 =  2.2;
// Inner radius of second straight section
      Float_t rInFaSaa1S2 = 17.9/2.;
// Length of second straight section
      Float_t dzFaSaa1S2  = 10.1;
// Length of flange SAA1 side      
//    Float_t dzFaSaa1F2  =  4.0;
// Inner radius of flange FA side
      Float_t rInFaSaa1F2 = 25.2/2.;
// Length of joint
      Float_t dzFaSaa1    = 34.8;
// Outer Radius at the end of the joint
      Float_t rOuFaSaa1E  = 41.93/2.;
      
      
      TGeoPcon* shFaSaa1 = new TGeoPcon(0., 360., 8);
      z = 0;
// Flange FA side      
      shFaSaa1->DefineSection( 0, z, rInFaSaa1F1, rOuFaSaa1F1);
      z += dzFaSaa1F1;
      shFaSaa1->DefineSection( 1, z, rInFaSaa1F1, 40.0);
      shFaSaa1->DefineSection( 2, z, rInFaSaa1S1, 40.0);
// First straight section
      z += dzFaSaa1S1;
      shFaSaa1->DefineSection( 3, z, rInFaSaa1S1, 40.0);
// 45 deg transition region
      z += dzFaSaa1T1;
      shFaSaa1->DefineSection( 4, z, rInFaSaa1S2, 40.0);
// Second straight section
      z += dzFaSaa1S2;
      shFaSaa1->DefineSection( 5, z, rInFaSaa1S2, 40.0);
      shFaSaa1->DefineSection( 6, z, rInFaSaa1F2, 40.0);
// Flange SAA1 side
      z = dzFaSaa1;
      shFaSaa1->DefineSection( 7, z, rInFaSaa1F2, rOuFaSaa1E);

// Outer 2 deg line
      for (Int_t i  = 1; i < 7; i++) {
	  Double_t z    = shFaSaa1->GetZ(i);
	  Double_t r1   = shFaSaa1->GetRmin(i);	 
	  Double_t r2 = 39.5/2. + z * TMath::Tan(2. * kDegRad) - 0.01;
	  shFaSaa1->DefineSection(i, z, r1, r2);
      }
      TGeoVolume* voFaSaa1 = new TGeoVolume("YFASAA1", shFaSaa1, kMedNiWsh);
//
// Outer region with lower transport cuts
      TGeoCone* shFaSaa1O   = new TGeoCone(dzFaSaa1/2., rOuFaSaa1F1 - 3.5, rOuFaSaa1F1, rOuFaSaa1E - 3.5, rOuFaSaa1E);
      TGeoVolume* voFaSaa1O = new TGeoVolume("YFASAA1O", shFaSaa1O, kMedNiW);
      voFaSaa1->AddNode(voFaSaa1O, 1, new TGeoTranslation(0., 0., dzFaSaa1/2.));
      

///////////////////////////////////
//    SAA1 Steel Envelope        //
//    Drawing ALIP2A__0039       //
///////////////////////////////////
      
      Float_t rOut;      // Outer radius
// Thickness of the steel envelope
      Float_t dSt = 4.;
// 4 Section
// z-positions
      Float_t zSaa1StEnv[5] = {111.2, 113.7, 229.3, 195.0};
// Radii
// 1
      Float_t rOuSaa1StEnv1 = 40.4/2.;
      Float_t rInSaa1StEnv1 = rOuSaa1StEnv1 - dSt;
// 2  
      Float_t rInSaa1StEnv2 = 41.7/2.;
      Float_t rOuSaa1StEnv2 = rInSaa1StEnv2 + dSt / TMath::Cos(2.0 * kDegRad);
// 3         
      Float_t rOuSaa1StEnv3 = 57.6/2.;
      Float_t rInSaa1StEnv3 = rOuSaa1StEnv3 - dSt;
// 4
      Float_t rInSaa1StEnv4 = 63.4/2.;
      Float_t rOuSaa1StEnv4 = rInSaa1StEnv4 + dSt / TMath::Cos(1.6 * kDegRad);
// end         
      Float_t rInSaa1StEnv5 = 74.28/2.;
      Float_t rOuSaa1StEnv5 = rInSaa1StEnv5 + dSt / TMath::Cos(1.6 * kDegRad);
// Relative starting position
      Float_t zSaa1StEnvS = 3.;
      
      TGeoPcon* shSaa1StEnv  = new TGeoPcon(0., 360., 11);
// 1st Section
      z = zSaa1StEnvS;
      shSaa1StEnv->DefineSection( 0, z, rInSaa1StEnv1, rOuSaa1StEnv1);
      z += (zSaa1StEnv[0] - dSt);
      shSaa1StEnv->DefineSection( 1, z, rInSaa1StEnv1, rOuSaa1StEnv1);
// 1 - 2
      shSaa1StEnv->DefineSection( 2, z, rInSaa1StEnv1, rOuSaa1StEnv2);
      z +=  dSt;
      shSaa1StEnv->DefineSection( 3, z, rInSaa1StEnv1, rOuSaa1StEnv2);      
// 2nd Section
      shSaa1StEnv->DefineSection( 4, z, rInSaa1StEnv2, rOuSaa1StEnv2);      
      z += zSaa1StEnv[1];
      shSaa1StEnv->DefineSection( 5, z, rInSaa1StEnv3, rOuSaa1StEnv3);      
// 3rd Section  
      z += (zSaa1StEnv[2] - dSt);
      shSaa1StEnv->DefineSection( 6, z, rInSaa1StEnv3, rOuSaa1StEnv3);      
// 3 - 4
      shSaa1StEnv->DefineSection( 7, z, rInSaa1StEnv3, rOuSaa1StEnv4);
      z += dSt;
      shSaa1StEnv->DefineSection( 8, z, rInSaa1StEnv3, rOuSaa1StEnv4);
// 4th Section
      shSaa1StEnv->DefineSection( 9, z, rInSaa1StEnv4, rOuSaa1StEnv4);
      z += zSaa1StEnv[3];
      shSaa1StEnv->DefineSection(10, z, rInSaa1StEnv5, rOuSaa1StEnv5);      
      TGeoVolume* voSaa1StEnv  =  new TGeoVolume("YSAA1_SteelEnvelope", shSaa1StEnv, kMedSteel);


      
 

///////////////////////////////////
//    SAA1 W-Pipe                //
//    Drawing ALIP2A__0059       //
///////////////////////////////////
//
//    Flange FA side
//    Length of first section      
      Float_t dzSaa1WPipeF1  = 0.9;
//    Outer radius 
      Float_t rOuSaa1WPipeF1 = 24.5/2.;
//    Inner Radius
      Float_t rInSaa1WPipeF1 = 22.0/2.;
//    Length of second section
      Float_t dzSaa1WPipeF11 =  2.1;
//    Inner Radius
      Float_t rInSaa1WPipeF11 = 18.5/2.;
//
//    Central tube 
//    Length 
      Float_t dzSaa1WPipeC = 111.2;
//    Inner Radius at the end
      Float_t rInSaa1WPipeC = 22.0/2.;
//    Outer Radius
      Float_t rOuSaa1WPipeC = 31.9/2.;
//
//    Flange SAA2 Side
//    Length
      Float_t dzSaa1WPipeF2  = 6.0;
//    Outer radius 
      Float_t rOuSaa1WPipeF2 = 41.56/2.;

//
      TGeoPcon* shSaa1WPipe  = new TGeoPcon(0., 360., 8);
      z = 0.;
// Flange FA side first section
      shSaa1WPipe->DefineSection( 0, z, rInSaa1WPipeF1 , rOuSaa1WPipeF1);
      z += dzSaa1WPipeF1;
      shSaa1WPipe->DefineSection( 1, z, rInSaa1WPipeF1 , rOuSaa1WPipeF1);
// Flange FA side second section
      shSaa1WPipe->DefineSection( 2, z, rInSaa1WPipeF11, rOuSaa1WPipeF1);
      z += dzSaa1WPipeF11;
      shSaa1WPipe->DefineSection( 3, z, rInSaa1WPipeF11, rOuSaa1WPipeF1);
// Central Section
      shSaa1WPipe->DefineSection( 4, z, rInSaa1WPipeF11, rOuSaa1WPipeC);
      z += dzSaa1WPipeC;
      shSaa1WPipe->DefineSection( 5, z, rInSaa1WPipeC,   rOuSaa1WPipeC);
// Flange SAA2 side
      shSaa1WPipe->DefineSection( 6, z, rInSaa1WPipeC,   rOuSaa1WPipeF2);
      z +=  dzSaa1WPipeF2;
      shSaa1WPipe->DefineSection( 7, z, rInSaa1WPipeC,   rOuSaa1WPipeF2);

      TGeoVolume* voSaa1WPipe  =  new TGeoVolume("YSAA1_WPipe", shSaa1WPipe, kMedNiW);
//
// Inner region with higher transport cuts
      TGeoTube*   shSaa1WPipeI = new TGeoTube(rInSaa1WPipeC, rOuSaa1WPipeC, dzSaa1WPipeC/2.);
      TGeoVolume* voSaa1WPipeI =  new TGeoVolume("YSAA1_WPipeI", shSaa1WPipeI, kMedNiWsh);
      voSaa1WPipe->AddNode(voSaa1WPipeI, 1, new TGeoTranslation(0., 0., dzSaa1WPipeF1 + dzSaa1WPipeF11 + dzSaa1WPipeC/2));
      

///////////////////////////////////
//    SAA1 Pb Components         //
//    Drawing ALIP2A__0078       //
///////////////////////////////////
//
//    Inner angle
      Float_t tanAlpha = TMath::Tan(1.69    / 2. * kDegRad);
      Float_t tanBeta  = TMath::Tan(3.20    / 2. * kDegRad);    
// 
//    1st Section 2deg opening cone
//    Length 
      Float_t dzSaa1PbComp1 = 100.23;
//    Inner radius at entrance
      Float_t rInSaa1PbComp1 = 22.0/2.; // It's 21 cm diameter in the drawing. Is this a typo ??!!
//    Outer radius at entrance
      Float_t rOuSaa1PbComp1 = 42.0/2.;
//
//    2nd Section: Straight Section
//    Length 
      Float_t dzSaa1PbComp2  = 236.77;
//    Inner radius
      Float_t rInSaa1PbComp2 = rInSaa1PbComp1 + dzSaa1PbComp1 * tanAlpha;
//    Outer radius
      Float_t rOuSaa1PbComp2 = 49.0/2.;
//
//    3rd Section: 1.6deg opening cone until bellow
//    Length 
      Float_t dzSaa1PbComp3 = 175.6;
//    Inner radius        
      Float_t  rInSaa1PbComp3 = rInSaa1PbComp2 + dzSaa1PbComp2 * tanAlpha;
//    Outer radius
      Float_t rOuSaa1PbComp3 = 62.8/2.;
//   
//   4th Section: Bellow region
      Float_t dzSaa1PbComp4  = 26.4;
//    Inner radius        
      Float_t  rInSaa1PbComp4 = 37.1/2.;
      Float_t  rInSaa1PbCompB = 43.0/2.;
//    Outer radius
      Float_t rOuSaa1PbComp4  = rOuSaa1PbComp3 +  dzSaa1PbComp3 * tanBeta;
//
//   5th Section: Flange SAA2 side
//   1st detail 
      Float_t dzSaa1PbCompF1  = 4.;
      Float_t rOuSaa1PbCompF1 = 74.1/2.;
//   2nd detail       
      Float_t dzSaa1PbCompF2  = 3.;
      Float_t rOuSaa1PbCompF2 = 66.0/2.;
      Float_t rOuSaa1PbCompF3 = 58.0/2.;
      

      TGeoPcon* shSaa1PbComp  = new TGeoPcon(0., 360., 11);
      z = 120.2;
// 2 deg opening cone
      shSaa1PbComp->DefineSection( 0, z, rInSaa1PbComp1, rOuSaa1PbComp1);
      z += dzSaa1PbComp1;
      shSaa1PbComp->DefineSection( 1, z, rInSaa1PbComp2, rOuSaa1PbComp2);
// Straight section
      z += dzSaa1PbComp2;
      shSaa1PbComp->DefineSection( 2, z, rInSaa1PbComp3, rOuSaa1PbComp2);
// 1.6 deg opening cone
      shSaa1PbComp->DefineSection( 3, z, rInSaa1PbComp3, rOuSaa1PbComp3);
      z += dzSaa1PbComp3;
      shSaa1PbComp->DefineSection( 4, z, rInSaa1PbComp4, rOuSaa1PbComp4);
// Bellow region until outer flange
      shSaa1PbComp->DefineSection( 5, z, rInSaa1PbCompB, rOuSaa1PbComp4);
      z += (dzSaa1PbComp4 - dzSaa1PbCompF1 -  dzSaa1PbCompF2);
      shSaa1PbComp->DefineSection( 6, z, rInSaa1PbCompB, rOuSaa1PbCompF1);
      shSaa1PbComp->DefineSection( 7, z, rInSaa1PbCompB, rOuSaa1PbCompF2);
// Flange first step
      z +=  dzSaa1PbCompF1;
      shSaa1PbComp->DefineSection( 8, z, rInSaa1PbCompB, rOuSaa1PbCompF2);
      shSaa1PbComp->DefineSection( 9, z, rInSaa1PbCompB, rOuSaa1PbCompF3);
// Flange second step
      z +=  dzSaa1PbCompF2;
      shSaa1PbComp->DefineSection( 10, z, rInSaa1PbCompB, rOuSaa1PbCompF3);

      TGeoVolume* voSaa1PbComp  =  new TGeoVolume("YSAA1_PbComp", shSaa1PbComp, kMedPb);
//
// Inner region with higher transport cuts
      TGeoPcon*   shSaa1PbCompI = MakeShapeFromTemplate(shSaa1PbComp, 0., -3.);
      TGeoVolume* voSaa1PbCompI =  new TGeoVolume("YSAA1_PbCompI", shSaa1PbCompI, kMedPbSh);
      voSaa1PbComp->AddNode(voSaa1PbCompI, 1, vec0);
      
///////////////////////////////////
//    SAA1 W-Cone                //
//    Drawing ALIP2A__0058       //
///////////////////////////////////
      // Length of the Cone
      Float_t dzSaa1WCone = 52.9;
      // Inner and outer radii
      Float_t rInSaa1WCone1 = 20.4;
      Float_t rOuSaa1WCone1 = rInSaa1WCone1 + 0.97;
      Float_t rOuSaa1WCone2 = rInSaa1WCone1 + 2.80;
      // relative z-position 
      Float_t zSaa1WCone    = 9.3;
      
      
      TGeoPcon* shSaa1WCone  = new TGeoPcon(0., 360., 2);
      z = zSaa1WCone;
      shSaa1WCone->DefineSection( 0, z, rInSaa1WCone1, rOuSaa1WCone1);
      z += dzSaa1WCone;
      shSaa1WCone->DefineSection( 1, z, rInSaa1WCone1, rOuSaa1WCone2);
      TGeoVolume* voSaa1WCone  =  new TGeoVolume("YSAA1_WCone", shSaa1WCone, kMedNiW);

///////////////////////////////////
//    SAA1 Steel-Ring            //
//    Drawing ALIP2A__0040       //
///////////////////////////////////
//
//    Length of the ring
      Float_t dzSaa1StRing = 4.;
//    Inner and outer radius
      Float_t rInSaa1String = 33.0;
      Float_t rOuSaa1String = 41.1;      
//    Relative z-position
      Float_t zSaa1StRing   = 652.2;
      TGeoPcon* shSaa1StRing  = new TGeoPcon(0., 360., 2);
      z = zSaa1StRing;
      shSaa1StRing->DefineSection( 0, z, rInSaa1String, rOuSaa1String);
      z += dzSaa1StRing;
      shSaa1StRing->DefineSection( 1, z, rInSaa1String, rOuSaa1String);
      TGeoVolume* voSaa1StRing  =  new TGeoVolume("YSAA1_StRing", shSaa1StRing, kMedSteel);

///////////////////////////////////
//    SAA1 Inner Tube            //
//    Drawing ALIP2A__0082       //
///////////////////////////////////
//
// Length of saa2:               659.2 cm
// Length of inner tube:         631.9 cm
// Lenth of bellow cavern:        27.3 cm    
// Radius at entrance 18.5/2,  d = 0.3
// Radius at exit     37.1/2,  d = 0.3 
//
      Float_t dzSaa1InnerTube     = 631.9/2.;  // Half length of the tube  
      Float_t rInSaa1InnerTube    =  18.2/2.;  // Radius at entrance
      Float_t rOuSaa1InnerTube    =  36.8/2.;  // Radius at exit
      Float_t dSaa1InnerTube      =   0.2   ;  // Thickness      
      TGeoVolume* voSaa1InnerTube = new TGeoVolume("YSAA1_InnerTube", 
						   new TGeoCone(dzSaa1InnerTube, 
								rInSaa1InnerTube - dSaa1InnerTube, rInSaa1InnerTube,
								rOuSaa1InnerTube - dSaa1InnerTube, rOuSaa1InnerTube),
						   kMedSteelSh);

///////////////////////////////////
//    SAA1 Outer Shape           //
//    Drawing ALIP2A__0107       //
///////////////////////////////////
      // Total length 
      Float_t dzSaa1 = 659.2;
      //
      TGeoPcon* shSaa1M  = new TGeoPcon(0., 360., 20);
      Float_t kSec = 0.01; // security distance to avoid trivial extrusions
      Float_t rmin = rInSaa1InnerTube - dSaa1InnerTube - kSec;
      rmax = rOuSaa1InnerTube - dSaa1InnerTube - kSec;
      z = 0.;
      shSaa1M->DefineSection( 0, z, rmin, rOuSaa1WPipeF1);
      z += dzSaa1WPipeF1;
      shSaa1M->DefineSection( 1, z, rmin, rOuSaa1WPipeF1);
      shSaa1M->DefineSection( 2, z, 0.,   rOuSaa1WPipeF1);
      z += dzSaa1WPipeF11;
      shSaa1M->DefineSection( 3, z, 0.,   rOuSaa1WPipeF1);
      shSaa1M->DefineSection( 4, z, 0.,   rOuSaa1StEnv1);
      z = zSaa1WCone;
      shSaa1M->DefineSection( 5, z, 0.,   rOuSaa1StEnv1);
      shSaa1M->DefineSection( 6, z, 0.,   rOuSaa1WCone1);
      z += dzSaa1WCone;
      shSaa1M->DefineSection( 7, z, 0.,   rOuSaa1WCone2);
      shSaa1M->DefineSection( 8, z, 0.,   rOuSaa1StEnv1);
      z =  zSaa1StEnv[0] - dSt + zSaa1StEnvS;
      shSaa1M->DefineSection( 9, z, 0.,   rOuSaa1StEnv1);
      shSaa1M->DefineSection(10, z, 0.,   rOuSaa1StEnv2);
      z +=  (zSaa1StEnv[1] + dSt);
      shSaa1M->DefineSection(11, z, 0.,   rOuSaa1StEnv3);
      z +=  (zSaa1StEnv[2] - dSt);
      shSaa1M->DefineSection(12, z, 0.,   rOuSaa1StEnv3);
      shSaa1M->DefineSection(13, z, 0.,   rOuSaa1StEnv4);

      z += (zSaa1StEnv[3] - dSt +  dzSaa1PbCompF1 + dzSaa1PbCompF2 - dzSaa1PbComp4);
      Float_t rmaxSaa1 = shSaa1M->GetRmax(13) + (z - shSaa1M->GetZ(13)) * TMath::Tan(1.6 * kDegRad);

      shSaa1M->DefineSection(14, z, 0.,    rmaxSaa1);
      shSaa1M->DefineSection(15, z, rmax,  rmaxSaa1);
      z = zSaa1StRing;     
      shSaa1M->DefineSection(16, z, rmax, rOuSaa1String);
      z += dzSaa1PbCompF1;
      shSaa1M->DefineSection(17, z, rmax, rOuSaa1String);
      shSaa1M->DefineSection(18, z, rmax, rOuSaa1PbCompF3);
      z += dzSaa1PbCompF2;
      shSaa1M->DefineSection(19, z, rmax, rOuSaa1PbCompF3);

//
//    Inner 1.69deg line
      for (Int_t i  = 2; i < 15; i++) {
	  Double_t z    = shSaa1M->GetZ(i);
	  Double_t r2 = shSaa1M->GetRmax(i);	 
	  Double_t r1 = rmin + (z - 0.9) * TMath::Tan(1.69 / 2. * kDegRad) - kSec;
	  shSaa1M->DefineSection(i, z, r1, r2);
      }

      TGeoVolume* voSaa1M  =  new TGeoVolume("YSAA1M", shSaa1M, kMedAir);


///////////////////////////////////
//                               //
// Recess Station 2              //
// Drawing ALIP2A__0260          //
///////////////////////////////////
///////////////////////////////////
//    SAA1 W-Ring 1              //
//    Drawing ALIP2A__0217       //
///////////////////////////////////
      Float_t saa1Wring1Width  =  5.85;
      TGeoPcon* shSaa1Wring1    = new TGeoPcon(0., 360., 2);
      shSaa1Wring1->DefineSection(0, 0.00           , 20.30, 23.175);
      shSaa1Wring1->DefineSection(1, saa1Wring1Width, 20.30, 23.400);
      TGeoVolume* voSaa1Wring1  =  new TGeoVolume("YSAA1_WRING1", shSaa1Wring1, kMedNiW);

///////////////////////////////////
//    SAA1 W-Ring 2              //
//    Drawing ALIP2A__0055       //
///////////////////////////////////
      Float_t saa1Wring2Rinner  = 20.30;
      Float_t saa1Wring2Router  = 23.40;
      Float_t saa1Wring2HWidth  =  3.75;
      Float_t saa1Wring2Cutoffx =  4.45;
      Float_t saa1Wring2Cutoffy =  4.45;
      TGeoTubeSeg* shSaa1Wring2a  = new TGeoTubeSeg(saa1Wring2Rinner, saa1Wring2Router, saa1Wring2HWidth, 0., 90.);
      shSaa1Wring2a->SetName("shSaa1Wring2a");
      TGeoBBox* shSaa1Wring2b  = new TGeoBBox(saa1Wring2Router / 2., saa1Wring2Router / 2., saa1Wring2HWidth);
      shSaa1Wring2b->SetName("shSaa1Wring2b");
      TGeoTranslation* trSaa1Wring2b 
	  = new TGeoTranslation("trSaa1Wring2b", saa1Wring2Router / 2. + saa1Wring2Cutoffx, saa1Wring2Router / 2. + saa1Wring2Cutoffy, 0.);
      trSaa1Wring2b->RegisterYourself();
      TGeoCompositeShape*  shSaa1Wring2 = new TGeoCompositeShape("shSaa1Wring2", "(shSaa1Wring2a)*(shSaa1Wring2b:trSaa1Wring2b)");
      TGeoVolume* voSaa1Wring2 = new TGeoVolume("YSAA1_WRING2", shSaa1Wring2, kMedNiW);

///////////////////////////////////
//    SAA1 W-Ring 3              //
//    Drawing ALIP2A__0216       //
///////////////////////////////////

      Float_t saa1Wring3Rinner  = 20.30;
      Float_t saa1Wring3Router  = 23.40;
      Float_t saa1Wring3HWidth  =  3.75;
      Float_t saa1Wring3Cutoffx =  4.50;
      Float_t saa1Wring3Cutoffy =  4.40;
      TGeoTubeSeg* shSaa1Wring3a  = new TGeoTubeSeg(saa1Wring3Rinner, saa1Wring3Router, saa1Wring3HWidth, 0., 90.);
      shSaa1Wring3a->SetName("shSaa1Wring3a");
      TGeoBBox* shSaa1Wring3b  = new TGeoBBox(saa1Wring3Router / 2., saa1Wring3Router / 2., saa1Wring3HWidth);
      shSaa1Wring3b->SetName("shSaa1Wring3b");
      TGeoTranslation* trSaa1Wring3b 
	  = new TGeoTranslation("trSaa1Wring3b", saa1Wring3Router / 2. + saa1Wring3Cutoffx, saa1Wring3Router / 2. + saa1Wring3Cutoffy, 0.);
      trSaa1Wring3b->RegisterYourself();
      TGeoCompositeShape*  shSaa1Wring3 = new TGeoCompositeShape("shSaa1Wring3", "(shSaa1Wring3a)*(shSaa1Wring3b:trSaa1Wring3b)");
      TGeoVolume* voSaa1Wring3 = new TGeoVolume("YSAA1_WRING3", shSaa1Wring3, kMedNiW);

///////////////////////////////////
//    SAA1 W-Ring 4              //
//    Drawing ALIP2A__0215       //
///////////////////////////////////
      Float_t saa1Wring4Width  =  5.85;
      TGeoPcon* shSaa1Wring4  = new TGeoPcon(0., 360., 5);
      shSaa1Wring4->DefineSection(0, 0.00, 20.30, 23.40);
      shSaa1Wring4->DefineSection(1, 1.00, 20.30, 23.40);
      shSaa1Wring4->DefineSection(2, 1.00, 20.30, 24.50);      
      shSaa1Wring4->DefineSection(3, 4.85, 20.30, 24.80);
      shSaa1Wring4->DefineSection(4, 5.85, 24.10, 24.80);
      TGeoVolume* voSaa1Wring4  =  new TGeoVolume("YSAA1_WRING4", shSaa1Wring4, kMedNiW);

///////////////////////////////////
//    SAA1 W-Ring 5              //
//    Drawing ALIP2A__0218       //
///////////////////////////////////
      Float_t saa1Wring5Rinner = 20.30;
      Float_t saa1Wring5Router = 23.40;
      Float_t saa1Wring5HWidth =  0.85;
      TGeoVolume* voSaa1Wring5    = new TGeoVolume("YSAA1_WRING5", 
						   new TGeoTube(saa1Wring5Rinner, saa1Wring5Router, saa1Wring5HWidth), kMedNiW);
//
// Position the rings in the assembly 
//      
      TGeoVolumeAssembly* asSaa1ExtraShield = new TGeoVolumeAssembly("YSAA1ExtraShield");
// Distance between rings
      Float_t saa1DWrings = 2.3;
//
      dz = - (saa1Wring1Width + 6. * saa1Wring2HWidth + 2. * saa1Wring3HWidth + saa1Wring4Width + 2. * saa1Wring5HWidth + 2. * saa1DWrings) / 2.;
      asSaa1ExtraShield->AddNode(voSaa1Wring1,    1, new TGeoTranslation(0., 0., dz));
      dz +=   saa1Wring1Width;
      dz +=   saa1Wring2HWidth;   
      asSaa1ExtraShield->AddNode(voSaa1Wring2,    1, new TGeoCombiTrans(0., 0., dz, rot000));
      asSaa1ExtraShield->AddNode(voSaa1Wring2,    2, new TGeoCombiTrans(0., 0., dz, rot180));
      dz +=   saa1Wring2HWidth;   
      dz +=   saa1DWrings;
      dz +=   saa1Wring2HWidth;   
      asSaa1ExtraShield->AddNode(voSaa1Wring2,    3, new TGeoCombiTrans(0., 0., dz, rot090));
      asSaa1ExtraShield->AddNode(voSaa1Wring2,    4, new TGeoCombiTrans(0., 0., dz, rot270));
      dz +=   saa1Wring2HWidth;   
      dz +=   saa1Wring5HWidth;   
      asSaa1ExtraShield->AddNode(voSaa1Wring5,    1, new TGeoTranslation(0., 0., dz));
      dz +=   saa1Wring5HWidth;   
      dz +=   saa1Wring2HWidth;   
      asSaa1ExtraShield->AddNode(voSaa1Wring2,    5, new TGeoCombiTrans(0., 0., dz, rot000));
      asSaa1ExtraShield->AddNode(voSaa1Wring2,    6, new TGeoCombiTrans(0., 0., dz, rot180));
      dz +=   saa1Wring2HWidth;   
      dz +=   saa1DWrings;
      dz +=   saa1Wring3HWidth;   
      asSaa1ExtraShield->AddNode(voSaa1Wring3,    1, new TGeoCombiTrans(0., 0., dz, rot090));
      asSaa1ExtraShield->AddNode(voSaa1Wring3,    2, new TGeoCombiTrans(0., 0., dz, rot270));
      dz +=   saa1Wring3HWidth;   
      asSaa1ExtraShield->AddNode(voSaa1Wring4,    1, new TGeoTranslation(0., 0., dz));
      dz +=   saa1Wring4Width;   
      const Float_t saa1ExtraShieldL = 48;
//
// Assemble SAA1
      voSaa1M->AddNode(voSaa1StEnv,     1, vec0);
      voSaa1M->AddNode(voSaa1WPipe,     1, vec0);
      voSaa1M->AddNode(voSaa1PbComp,    1, vec0);
      voSaa1M->AddNode(voSaa1WCone,     1, vec0);
      voSaa1M->AddNode(voSaa1StRing,    1, vec0);
      voSaa1M->AddNode(voSaa1InnerTube, 1, new TGeoTranslation(0., 0., dzSaa1InnerTube + 0.9));		      
      TGeoVolumeAssembly* voSaa1 = new TGeoVolumeAssembly("YSAA1");
      voSaa1->AddNode(voSaa1M, 1, vec0);
      
///////////////////////////////////////
//          SAA1/SAA2  Pb Joint      //
//          Drawing ALIP2A__0081     //
///////////////////////////////////////
//
// Outer radius
      Float_t rOuSaa1Saa2      = 70.0/2.;
// Flange SAA1 side
      Float_t dzSaa1Saa2F1     =  3.;
      Float_t rInSaa1Saa2F1    = 58.5/2.;
// 1st Central Section
      Float_t dzSaa1Saa2C1     = 19.3;
      Float_t rInSaa1Saa2C1    = 42.8/2.;
// Transition Region
      Float_t dzSaa1Saa2T      =  3.3;
// 1st Central Section
      Float_t dzSaa1Saa2C2     =  6.2;
      Float_t rInSaa1Saa2C2    = 36.2/2.;
// Flange SAA2 side
      Float_t dzSaa1Saa2F2     =  3.1;
      Float_t rInSaa1Saa2F2    = 54.1/2.;
// Total length
      Float_t dzSaa1Saa2       = 34.9;
      
      
      TGeoPcon* shSaa1Saa2Pb = new TGeoPcon(0., 360., 8);
      z = 0.;
// Flange SAA1 side
      shSaa1Saa2Pb->DefineSection( 0,  z, rInSaa1Saa2F1, rOuSaa1Saa2);
      z += dzSaa1Saa2F1;
      shSaa1Saa2Pb->DefineSection( 1,  z, rInSaa1Saa2F1, rOuSaa1Saa2);
      shSaa1Saa2Pb->DefineSection( 2,  z, rInSaa1Saa2C1, rOuSaa1Saa2);
// Central region 1
      z +=  dzSaa1Saa2C1;
      shSaa1Saa2Pb->DefineSection( 3,  z, rInSaa1Saa2C1, rOuSaa1Saa2);
// 45 deg transition
      z += dzSaa1Saa2T;
      shSaa1Saa2Pb->DefineSection( 4, z, rInSaa1Saa2C2, rOuSaa1Saa2);
      z += dzSaa1Saa2C2;
      shSaa1Saa2Pb->DefineSection( 5, z, rInSaa1Saa2C2, rOuSaa1Saa2);
      shSaa1Saa2Pb->DefineSection( 6, z, rInSaa1Saa2F2, rOuSaa1Saa2);
      z += dzSaa1Saa2F2;
      shSaa1Saa2Pb->DefineSection( 7, z, rInSaa1Saa2F2, rOuSaa1Saa2);
      TGeoVolume* voSaa1Saa2Pb = new TGeoVolume("YSAA1SAA2Pb", shSaa1Saa2Pb, kMedPb);
//
//    Mother volume and outer steel envelope
      Float_t rOuSaa1Saa2Steel = 36.9;
      
      TGeoPcon* shSaa1Saa2 = MakeShapeFromTemplate(shSaa1Saa2Pb, 0., rOuSaa1Saa2Steel-rOuSaa1Saa2);
      TGeoVolume* voSaa1Saa2 = new TGeoVolume("YSAA1SAA2", shSaa1Saa2, kMedSteel);
      voSaa1Saa2->AddNode(voSaa1Saa2Pb, 1, vec0);
//
//    Inner region with higher transport cuts
//
      TGeoPcon*   shSaa1Saa2I = MakeShapeFromTemplate(shSaa1Saa2Pb, 0., -3.);
      TGeoVolume* voSaa1Saa2I = new TGeoVolume("YSAA1_SAA2I", shSaa1Saa2I, kMedPbSh);
      voSaa1Saa2Pb->AddNode(voSaa1Saa2I, 1, vec0);
      


///////////////////////////////////////
//                SAA2               //
///////////////////////////////////////

      
///////////////////////////////////
//    SAA2 Steel Envelope        //
//    Drawing ALIP2A__0041       //
///////////////////////////////////
      dSt = 4.;  // Thickness of steel envelope
      // Length of the first section
      Float_t dzSaa2StEnv1  = 163.15;
      Float_t rInSaa2StEnv1 = 65.8/2.;
      // Length of the second section
      Float_t dzSaa2StEnv2  = 340.35;
      Float_t rInSaa2StEnv2 = 87.2/2.;
      // Rel. starting position 
      Float_t zSaa2StEnv = 3.;
      
      TGeoPcon* shSaa2StEnv  = new TGeoPcon(0., 360., 6);
      // First Section
      z = zSaa2StEnv;
      shSaa2StEnv->DefineSection( 0, z, rInSaa2StEnv1, rInSaa2StEnv1 + dSt);
      z += dzSaa2StEnv1;
      shSaa2StEnv->DefineSection( 1, z, rInSaa2StEnv1, rInSaa2StEnv1 + dSt);
      // Transition region
      shSaa2StEnv->DefineSection( 2, z, rInSaa2StEnv1, rInSaa2StEnv2 + dSt);           
      z += dSt;
      shSaa2StEnv->DefineSection( 3, z, rInSaa2StEnv1, rInSaa2StEnv2 + dSt);          
      // Second section
      shSaa2StEnv->DefineSection( 4, z, rInSaa2StEnv2, rInSaa2StEnv2 + dSt);           
      z += dzSaa2StEnv2;
      shSaa2StEnv->DefineSection( 5, z, rInSaa2StEnv2, rInSaa2StEnv2 + dSt);           

      TGeoVolume* voSaa2StEnv  =  new TGeoVolume("YSAA2_SteelEnvelope", shSaa2StEnv, kMedSteel);


///////////////////////////////////
//    SAA2 Pb Ring               //
//    Drawing ALIP2A__0080       //
//    Drawing ALIP2A__0111       //
///////////////////////////////////
//
// Rel. position in z
      Float_t zSaa2PbRing    = 35.25;  
// Length
      Float_t dzSaa2PbRing   = 65.90;  
// Inner radius
      Float_t rInSaa2PbRing  = 37.00;
// Outer radius at front
      Float_t rOuSaa2PbRingF = 42.74;
// Outer Rradius at rear
      Float_t rOuSaa2PbRingR = 44.58;      
      
      TGeoPcon* shSaa2PbRing  = new TGeoPcon(0., 360., 2);
      z = zSaa2PbRing;
      shSaa2PbRing->DefineSection(0, z, rInSaa2PbRing, rOuSaa2PbRingF);
      z += dzSaa2PbRing;
      shSaa2PbRing->DefineSection(1, z, rInSaa2PbRing, rOuSaa2PbRingR);

      TGeoVolume* voSaa2PbRing  =  new TGeoVolume("YSAA2_PbRing", shSaa2PbRing, kMedPb);


///////////////////////////////////
//    SAA2 Pb Components         //
//    Drawing ALIP2A__0079       //
///////////////////////////////////
      tanAlpha = TMath::Tan(1.89 / 2. * kDegRad);
      TGeoPcon* shSaa2PbComp  = new TGeoPcon(0., 360., 16);
      // Total length 
      Float_t dzSaa2PbComp    = 512.;
      // Length of 1st bellow recess
      Float_t dzSaa2PbCompB1  =  24.;
      // Length of 2nd bellow recess
      Float_t dzSaa2PbCompB2   = 27.;
      // Flange on the SAA1 side Detail A
      // 1st Step
      Float_t dzSaa2PbCompA1  =   1.5;
      Float_t rInSaa2PbCompA1 =  43.0/2.;
      Float_t rOuSaa2PbCompA1 =  53.0/2.;
      // 2nd Step
      Float_t dzSaa2PbCompA2  =   1.5;
      Float_t rInSaa2PbCompA2 =  36.8/2.;
      Float_t rOuSaa2PbCompA2 =  rOuSaa2PbCompA1;
      // Straight section
      Float_t dzSaa2PbCompA3  =  21.0;
      Float_t rInSaa2PbCompA3 =  rInSaa2PbCompA2;
      Float_t rOuSaa2PbCompA3 =  65.2/2.;
      //
      // 1st Section (outer straight, inner 1.89/2. deg opening cone)
      // Length
      Float_t dzSaa2PbComp1   = 146.15;
      // Inner radius at the end 
      Float_t rInSaa2PbComp1  = rInSaa2PbCompA3 + dzSaa2PbComp1 * tanAlpha;
      // Outer radius
      Float_t rOuSaa2PbComp1  = rOuSaa2PbCompA3;
      //
      // 2nd Section (outer straight, inner 1.89/2. deg opening cone)
      // Length 
      Float_t dzSaa2PbComp2   = (dzSaa2PbComp - dzSaa2PbComp1 - dzSaa2PbCompB1 - dzSaa2PbCompB2);
      // Inner radius at the end 
      Float_t rInSaa2PbComp2  = rInSaa2PbComp1 + dzSaa2PbComp2 * tanAlpha;
      // Outer radius
      Float_t rOuSaa2PbComp2  = 86.6/2.;
      //
      // Flange on the SAA3 side (Detail E)
      //
      // Straight Section
      // Length  dzSaa2PbCompB2 - 8.8 = 27 - 8.8 = 18.2
      Float_t dzSaa2PbCompE1  =  18.2;
      Float_t rInSaa2PbCompE1 =  52.0/2.;
      Float_t rOuSaa2PbCompE1 =  86.6/2.;
      // 45 deg transition
      Float_t dzSaa2PbCompE2  =   2.7;
      // 1st Step
      Float_t dzSaa2PbCompE3  =   0.6;
      Float_t rInSaa2PbCompE3 =  52.0/2.+ dzSaa2PbCompE2;
      Float_t rOuSaa2PbCompE3 =  83.0/2.;
      // 2nd Step
      Float_t dzSaa2PbCompE4  =   4.0;
      Float_t rOuSaa2PbCompE4 =  61.6/2.;
      // end
      Float_t dzSaa2PbCompE5  =   1.5;


      //
      // Flange on SAA1 side (Detail A)
      z = 0.;
      // 1st Step
      shSaa2PbComp->DefineSection(  0,  z, rInSaa2PbCompA1, rOuSaa2PbCompA1);
      z += dzSaa2PbCompA1;
      shSaa2PbComp->DefineSection(  1,  z, rInSaa2PbCompA1, rOuSaa2PbCompA1);
      shSaa2PbComp->DefineSection(  2,  z, rInSaa2PbCompA2, rOuSaa2PbCompA2);
      // 2nd Step
      z += dzSaa2PbCompA2;
      shSaa2PbComp->DefineSection(  3,  z, rInSaa2PbCompA2, rOuSaa2PbCompA2);
      shSaa2PbComp->DefineSection(  4,  z, rInSaa2PbCompA3, rOuSaa2PbCompA3);
      // straight section
      z += dzSaa2PbCompA3;
      shSaa2PbComp->DefineSection(  5,  z, rInSaa2PbCompA3, rOuSaa2PbCompA3);
      //
      // Section 1
      z += dzSaa2PbComp1;
      shSaa2PbComp->DefineSection(  6,  z, rInSaa2PbComp1, rOuSaa2PbComp1);
      //
      // Section 2
      shSaa2PbComp->DefineSection(  7,  z, rInSaa2PbComp2, rOuSaa2PbComp2);
      z += dzSaa2PbComp2;
      shSaa2PbComp->DefineSection(  8,  z, rInSaa2PbComp2, rOuSaa2PbComp2);
      //
      // Flange SAA3 side (Detail E)
      z += dzSaa2PbCompE1;
      shSaa2PbComp->DefineSection(  9,  z, rInSaa2PbCompE1, rOuSaa2PbCompE1);
      // 45 deg transition
      z += dzSaa2PbCompE2;
      shSaa2PbComp->DefineSection( 10,  z, rInSaa2PbCompE3, rOuSaa2PbCompE1);
      // 1st step
      z += dzSaa2PbCompE3;
      shSaa2PbComp->DefineSection( 11,  z, rInSaa2PbCompE3, rOuSaa2PbCompE1);
      shSaa2PbComp->DefineSection( 12,  z, rInSaa2PbCompE3, rOuSaa2PbCompE3);
      // 2nd step
      z += dzSaa2PbCompE4;
      shSaa2PbComp->DefineSection( 13,  z, rInSaa2PbCompE3, rOuSaa2PbCompE3);
      shSaa2PbComp->DefineSection( 14,  z, rInSaa2PbCompE3, rOuSaa2PbCompE4);
      // end
      z += dzSaa2PbCompE5;
      shSaa2PbComp->DefineSection( 15,  z, rInSaa2PbCompE3, rOuSaa2PbCompE4);

      TGeoVolume* voSaa2PbComp  =  new TGeoVolume("YSAA2_PbComp", shSaa2PbComp, kMedPb);


///////////////////////////////////
//    SAA2 Inner Tube            //
//    Drawing ALIP2A__0083       //
///////////////////////////////////
//
//
//
// Length of saa2:               512.0 cm
// Length of inner tube:         501.7 cm 
// Lenth of bellow recess:        10.3 cm   ( 1.5 + 8.8) 
// Radius at entrance 36.8/2,  d = 0.1
// Radius at exit     52.0/2,  d = 0.1 
//
      const Float_t kSaa2InnerTubeL     =  501.7;    // Length of the tube  
      const Float_t kSaa2InnerTubeRmin  =  36.6/2.;  // Radius at entrance
      const Float_t kSaa2InnerTubeRmax  =  51.8/2.;  // Radius at exit
      const Float_t kSaa2InnerTubeD     =   0.2   ;  // Thickness     
      TGeoPcon*   shSaa2InnerTube = new TGeoPcon(0., 360., 4);
      z = 0.;
      shSaa2InnerTube->DefineSection( 0, z, kSaa2InnerTubeRmin - kSaa2InnerTubeD, kSaa2InnerTubeRmin);
      z += dzSaa2PbCompA2 + dzSaa2PbCompA3;
      shSaa2InnerTube->DefineSection( 1, z, kSaa2InnerTubeRmin - kSaa2InnerTubeD, kSaa2InnerTubeRmin);
      z =  kSaa2InnerTubeL - dzSaa2PbCompE1;
      shSaa2InnerTube->DefineSection( 2, z, kSaa2InnerTubeRmax - kSaa2InnerTubeD, kSaa2InnerTubeRmax);
      z =  kSaa2InnerTubeL;
      shSaa2InnerTube->DefineSection( 3, z, kSaa2InnerTubeRmax - kSaa2InnerTubeD, kSaa2InnerTubeRmax);
      TGeoVolume* voSaa2InnerTube = new TGeoVolume("YSAA2_InnerTube", shSaa2InnerTube, kMedSteel);
      
///////////////////////////////////
//    SAA2 Steel Ring            //
//    Drawing ALIP2A__0042       //
///////////////////////////////////
      //  HalfWidth
      Float_t dzSaa2SteelRing = 2.;
      TGeoTube*   shSaa2SteelRing  = new TGeoTube(41.6, 47.6, dzSaa2SteelRing);
      TGeoVolume* voSaa2SteelRing  = new TGeoVolume("YSAA2_SteelRing", shSaa2SteelRing, kMedSteel);

///////////////////////////////////
//    SAA2 Outer Shape           //
//    Drawing ALIP2A__0108       //
///////////////////////////////////

      TGeoPcon* shSaa2  = new TGeoPcon(0., 360., 16);
      kSec = 0.02; // security distance to avoid trivial extrusions
      rmin = kSaa2InnerTubeRmin - kSaa2InnerTubeD - kSec;
      rmax = kSaa2InnerTubeRmax - kSaa2InnerTubeD - kSec;
      // Flange SAA1 side
      z = 0.;
      shSaa2->DefineSection( 0, z, rmin           , rOuSaa2PbCompA1);
      z += dzSaa2PbCompA1 + dzSaa2PbCompA2;
      shSaa2->DefineSection( 1, z, rmin           , rOuSaa2PbCompA1);
      shSaa2->DefineSection( 2, z, rmin           , rInSaa2StEnv1 + dSt);
      z += dzSaa2PbCompA3;
      shSaa2->DefineSection( 3, z, rmin           , rInSaa2StEnv1 + dSt);
      z = zSaa2PbRing;
      shSaa2->DefineSection( 4, z, 0.             , rInSaa2StEnv1 + dSt);
      shSaa2->DefineSection( 5, z, 0.             , rOuSaa2PbRingF);
      z += dzSaa2PbRing;
      shSaa2->DefineSection( 6, z, 0.             , rOuSaa2PbRingR);
      shSaa2->DefineSection( 7, z, 0.             , rInSaa2StEnv1 + dSt);
      z = dzSaa2PbCompA1 + dzSaa2PbCompA2 + dzSaa2StEnv1;
      shSaa2->DefineSection( 8, z, 0.             , rInSaa2StEnv1 + dSt);
      shSaa2->DefineSection( 9, z, 0.             , rInSaa2StEnv2 + dSt);
      z = dzSaa2PbComp - dzSaa2PbCompB2;
      shSaa2->DefineSection(10, z, rmax           , rInSaa2StEnv2 + dSt);
      z += dzSaa2PbCompE1;
      shSaa2->DefineSection(11, z, rmax           , rInSaa2StEnv2 + dSt);
      z += dzSaa2PbCompE2;
      shSaa2->DefineSection(12, z, rInSaa2PbCompE3, rInSaa2StEnv2 + dSt);
      z += (dzSaa2PbCompE3 + dzSaa2PbCompE4);
      shSaa2->DefineSection(13, z, rInSaa2PbCompE3, rInSaa2StEnv2 + dSt);
      shSaa2->DefineSection(14, z, rInSaa2PbCompE3, rOuSaa2PbCompE4);
      z += dzSaa2PbCompE5;
      shSaa2->DefineSection(15, z, rInSaa2PbCompE3, rOuSaa2PbCompE4);

      TGeoVolume* voSaa2  =  new TGeoVolume("YSAA2", shSaa2, kMedAir);
// Inner 1.89/2 deg line
      Double_t zref   = dzSaa2PbCompA1 + dzSaa2PbCompA2 + dzSaa2PbCompA3;
      for (Int_t i  = 4; i < 10; i++) {
	  Double_t z    = shSaa2->GetZ(i);
	  Double_t r2 = shSaa2->GetRmax(i);	 
	  Double_t r1 = rmin + (z -  zref) * TMath::Tan(1.89 / 2. * kDegRad) - kSec;
	  shSaa2->DefineSection(i, z, r1, r2);
      }

//
//    Assemble SAA2
      voSaa2->AddNode(voSaa2StEnv,     1, vec0);
      voSaa2->AddNode(voSaa2PbRing,    1, vec0);
      voSaa2->AddNode(voSaa2PbComp,    1, vec0);
      voSaa2->AddNode(voSaa2InnerTube, 1, new TGeoTranslation(0., 0., dzSaa2PbCompA1));
      z = (dzSaa2PbComp - dzSaa2PbCompE4 - dzSaa2PbCompE5) + dzSaa2SteelRing;
      voSaa2->AddNode(voSaa2SteelRing, 1, new TGeoTranslation(0., 0., z));


///////////////////////////////////////
//                SAA3               //
///////////////////////////////////////
//
//
//  This is a study performed by S. Maridor
//  The SAA3 has not yet been designed !!!!!!!!
//
///////////////////////////////////
//    SAA3 Outer Shape           //
//    Drawing ALIP2A__0xxx       //
///////////////////////////////////
      
      TGeoVolumeAssembly* voSaa3  = new TGeoVolumeAssembly("YSAA3");
      
///////////////////////////////////
//    SAA3 Steel Components      //
//    Drawing ALIP2A__0xxx       //
///////////////////////////////////
//    Block
      TGeoBBox* shSaa3SteelBlockO   = new TGeoBBox(220./2., 80./2., 100./2.);
      shSaa3SteelBlockO->SetName("Saa3SteelBlockO");

      TGeoPcon* shSaa3InnerRegion  = new TGeoPcon(0., 360., 6);
      shSaa3InnerRegion->DefineSection( 0, -52.0, 0., 56.6/2.);
      shSaa3InnerRegion->DefineSection( 1, -45.0, 0., 56.6/2.);
      shSaa3InnerRegion->DefineSection( 2, -42.0, 0., 50.6/2.);
      shSaa3InnerRegion->DefineSection( 3, -30.0, 0., 50.6/2.);
      shSaa3InnerRegion->DefineSection( 4,  30.5, 0., 16.8/2.);
      shSaa3InnerRegion->DefineSection( 5,  52.0, 0., 16.8/2.);
      shSaa3InnerRegion->SetName("Saa3InnerRegion");

      TGeoCompositeShape* shSaa3SteelBlock = new TGeoCompositeShape("Saa3SteelBlock", "Saa3SteelBlockO-Saa3InnerRegion");
      TGeoVolume* voSaa3SteelBlock         = new TGeoVolume("YSAA3SteelBlock", shSaa3SteelBlock, kMedSteel);     


      voSaa3->AddNode(voSaa3SteelBlock, 1, vec0);
//    Plate 1: 220 cm x 100 cm x 10 cm (x 1)
      TGeoBBox* shSaa3SteelPlate1   = new TGeoBBox(220./2., 10./2., 100./2.);
      TGeoVolume* voSaa3SteelPlate1  =  new TGeoVolume("YSAA3SteelPlate1", shSaa3SteelPlate1, kMedSteel);
      voSaa3->AddNode(voSaa3SteelPlate1, 1, new TGeoTranslation(0., -80./2. - 10./2., 0.));
//    Plate 2: 200 cm x 100 cm x 10 cm (x 6)
      TGeoBBox* shSaa3SteelPlate2   = new TGeoBBox(200./2., 60./2., 100./2.);
      TGeoVolume* voSaa3SteelPlate2  =  new TGeoVolume("YSAA3SteelPlate2", shSaa3SteelPlate2, kMedSteel);
      voSaa3->AddNode(voSaa3SteelPlate2, 1, new TGeoTranslation(0., -80./2. - 10. - 60./2., 0.));
//    Plate 3: 200 cm x 100 cm x 10 cm (x 4)
      TGeoBBox* shSaa3SteelPlate3   = new TGeoBBox(200./2., 40./2., 100./2.);
      TGeoVolume* voSaa3SteelPlate3  =  new TGeoVolume("YSAA3SteelPlate3", shSaa3SteelPlate3, kMedSteel);
      voSaa3->AddNode(voSaa3SteelPlate3, 1, new TGeoTranslation(0., +80./2. + 40/2., 0.));
//    Plate 4: 180 cm x 100 cm x 10 cm (x 1)
      TGeoBBox* shSaa3SteelPlate4   = new TGeoBBox(180./2., 10./2., 100./2.);
      TGeoVolume* voSaa3SteelPlate4  =  new TGeoVolume("YSAA3SteelPlate4", shSaa3SteelPlate4, kMedSteel);
      voSaa3->AddNode(voSaa3SteelPlate4, 1, new TGeoTranslation(0., +80./2. + 40. + 10/2., 0.));
//    Plate 5: 150 cm x 100 cm x 10 cm (x 1)
      TGeoBBox* shSaa3SteelPlate5   = new TGeoBBox(150./2., 10./2., 100./2.);
      TGeoVolume* voSaa3SteelPlate5  =  new TGeoVolume("YSAA3SteelPlate5", shSaa3SteelPlate5, kMedSteel);
      voSaa3->AddNode(voSaa3SteelPlate5, 1, new TGeoTranslation(0., +80./2. + 50. + 10/2., 0.));
//    Plate 6: 100 cm x 100 cm x 10 cm (x 1)
      TGeoBBox* shSaa3SteelPlate6   = new TGeoBBox(100./2., 10./2., 100./2.);
      TGeoVolume* voSaa3SteelPlate6  =  new TGeoVolume("YSAA3SteelPlate6", shSaa3SteelPlate6, kMedSteel);
      voSaa3->AddNode(voSaa3SteelPlate6, 1, new TGeoTranslation(0., +80./2. + 60. + 10/2., 0.));

///////////////////////////////////
//    SAA3 Concrete Piece        //
//    Drawing ALIP2A__0xxx       //
///////////////////////////////////
      TGeoPcon* shSaa3ConcPiece  = new TGeoPcon(0., 360., 6);
      shSaa3ConcPiece->DefineSection( 0,   0.0, 57.4/2., 62./2.);
      shSaa3ConcPiece->DefineSection( 1,   5.0, 57.4/2., 62./2.);
      shSaa3ConcPiece->DefineSection( 2,   8.0, 51.4/2., 62./2.);
      shSaa3ConcPiece->DefineSection( 3,  20.0, 51.4/2., 62./2.);
      shSaa3ConcPiece->DefineSection( 4,  80.5, 17.6/2., 62./2.);
      shSaa3ConcPiece->DefineSection( 5, 100.0, 17.6/2., 62./2.);
      TGeoVolume* voSaa3ConcPiece  =  new TGeoVolume("YSAA2_ConcPiece", shSaa3ConcPiece, kMedConc);
      voSaa3SteelBlock->AddNode(voSaa3ConcPiece, 1, new TGeoTranslation(0., 0., -50.));

///////////////////////////////////
//    SAA3 InnerTube             //
//    Drawing ALIP2A__0xxx       //
///////////////////////////////////
      TGeoPcon* shSaa3InnerTube  = new TGeoPcon(0., 360., 6);
      shSaa3InnerTube->DefineSection( 0,   0.0, 56.6/2., 57.0/2.);
      shSaa3InnerTube->DefineSection( 1,   5.0, 56.6/2., 57.0/2.);
      shSaa3InnerTube->DefineSection( 2,   8.0, 50.6/2., 51.0/2.);
      shSaa3InnerTube->DefineSection( 3,  20.0, 50.6/2., 51.0/2.);
      shSaa3InnerTube->DefineSection( 4,  80.5, 16.8/2., 17.2/2.);
      shSaa3InnerTube->DefineSection( 5, 100.0, 16.8/2., 17.2/2.);
      TGeoVolume* voSaa3InnerTube  =  new TGeoVolume("YSAA2_InnerTube", shSaa3InnerTube, kMedSteel);
      voSaa3SteelBlock->AddNode(voSaa3InnerTube, 1, new TGeoTranslation(0., 0., -50.));

///////////////////////////////////
//    Muon Filter                //
//    Drawing ALIP2A__0105       //
///////////////////////////////////
      // Half Length 
      Float_t dzMuonFilter = 60.;
      
      TGeoBBox*   shMuonFilterO  = new TGeoBBox(550./2., 620./2., dzMuonFilter);
      shMuonFilterO->SetName("FilterO");
      TGeoTube*   shMuonFilterI  = new TGeoTube(0., 50., dzMuonFilter);
      shMuonFilterI->SetName("FilterI");
      TGeoCompositeShape* shMuonFilter = new TGeoCompositeShape("MuonFilter", "FilterO-FilterI");
      //
      // !!!!! Needs to be inclined
      TGeoVolume* voMuonFilter = new TGeoVolume("YMuonFilter", shMuonFilter, kMedSteel);

//  
      TGeoVolumeAssembly* voSaa  = new TGeoVolumeAssembly("YSAA");
//
//    
//    
//
//
//    Starting position of the FA Flange/Tail 
      Float_t ziFaWTail   = 499.0;
//    End of the FA Flange/Tail
      Float_t zoFaWTail   = ziFaWTail + dzFaWTail;
//    Starting position of the FA/SAA1 Joint (2.8 cm overlap with tail)
      Float_t ozFaSaa1    = 2.8;
      Float_t ziFaSaa1    = zoFaWTail - ozFaSaa1;
//    End of the FA/SAA1 Joint
      Float_t zoFaSaa1    = ziFaSaa1 +  dzFaSaa1;
//    Starting position of SAA1 (2.0 cm overlap with joint)     
      Float_t ozSaa1      = 2.;
      Float_t ziSaa1      = zoFaSaa1 - ozSaa1;
//    End of SAA1
      Float_t zoSaa1      = ziSaa1 + dzSaa1;
//    Starting position of SAA1/SAA2 Joint (1.95 cm overlap with SAA1)
      Float_t ziSaa1Saa2  = zoSaa1 - 1.95;
//    End of SAA1/SAA2 Joint
      Float_t zoSaa1Saa2  = ziSaa1Saa2 + dzSaa1Saa2;
//    Starting position of SAA2 (3.1 cm overlap with the joint)
      Float_t ziSaa2      = zoSaa1Saa2 - 3.1;
//    End of SAA2
      Float_t zoSaa2      = ziSaa2 + dzSaa2PbComp;
//    Position of SAA3
      Float_t zcSaa3      = zoSaa2 + 50.;
//    Position of the Muon Filter
      Float_t zcFilter    = 1465.9 + dzMuonFilter;

      printf("FaSaa1    %f %f \n", ziFaSaa1, zoFaSaa1);
      printf("Saa1      %f %f \n", ziSaa1, zoSaa1);
      printf("Saa1/Saa2 %f %f \n", ziSaa1Saa2, zoSaa1Saa2);
      printf("Saa2      %f %f \n", ziSaa2, zoSaa2);
      
      voSaa->AddNode(voFaWTail,    1, new TGeoTranslation(0., 0., ziFaWTail));
      voSaa->AddNode(voFaSaa1,     1, new TGeoTranslation(0., 0., ziFaSaa1));
      voSaa->AddNode(voSaa1 ,      1, new TGeoTranslation(0., 0., ziSaa1));
      voSaa->AddNode(voSaa1Saa2,   1, new TGeoTranslation(0., 0., ziSaa1Saa2));
      voSaa->AddNode(voSaa2 ,      1, new TGeoTranslation(0., 0., ziSaa2));
      voSaa->AddNode(voSaa3,       1, new TGeoTranslation(0., 0., zcSaa3));
      voSaa->AddNode(voMuonFilter, 1, new TGeoTranslation(0., 0., zcFilter));      

      TGeoRotation* rotxz  = new TGeoRotation("rotxz",   90.,   0., 90.,  90., 180., 0.);
      top->AddNode(voSaa, 1, new TGeoCombiTrans(0., 0., 0., rotxz));
//
//  Mother volume for muon stations 1+2 and shielding material placed between the quadrants
//
      // Position of the dipole
      Float_t ziDipole = 724.45;
      
      TGeoPcon* shYOUT1 = new TGeoPcon(0., 360., 25);
      Float_t eps = 1.e-2;
      // FA Tail Section
      for (Int_t iz =  0; iz < 9; iz++) {
	  z = shFaWTail->GetZ(iz+1);
	  if (iz == 8) z -= ozFaSaa1;
	  shYOUT1->DefineSection(iz, z + ziFaWTail, shFaWTail->GetRmax(iz+1) + eps, 150.);
      }
      // FA-SAA1 Joint
      z = shYOUT1->GetZ(8);
      
      for (Int_t iz =  9; iz < 17; iz++) 
	  shYOUT1->DefineSection(iz, z + shFaSaa1->GetZ(iz-9), shFaSaa1->GetRmax(iz-9) + eps, 150.);

      z = shYOUT1->GetZ(16) - ozSaa1;
      // SAA1  - Dipole
      for (Int_t iz = 17; iz < 24; iz++) 
	  shYOUT1->DefineSection(iz, z + shSaa1M->GetZ(iz-13), shSaa1M->GetRmax(iz-13) + eps, 150.); 
      // Distance between dipole and start of SAA1 2deg opening cone
      dz   = zSaa1StEnv[0] - dSt + zSaa1StEnvS + ziSaa1 - ziDipole;
      rOut = rOuSaa1StEnv2 + dz * TMath::Tan(2. * kDegRad);
      
      shYOUT1->DefineSection(24, ziDipole, rOut + eps, 150.);

      InvertPcon(shYOUT1);
      TGeoVolume* voYOUT1 = new TGeoVolume("YOUT1", shYOUT1, kMedAirMu);

      voYOUT1->AddNode(asSaa1ExtraShield, 1, new TGeoCombiTrans(0., 0., - (100.7 + 62.2 + saa1ExtraShieldL / 2. + ziFaWTail), rotxz));
      voYOUT1->AddNode(asFaExtraShield,   1, new TGeoCombiTrans(0., 0., - (16.41 + kFaWring2HWidth + ziFaWTail), rotxz));
      top->AddNode(voYOUT1, 1, vec0);
//
//  Mother volume for muon stations 4+5 and trigger stations.
//
      Float_t zoDipole = 1235.55;
      
      TGeoPcon* shYOUT2 = new TGeoPcon(0., 360., 14);
      z =  zoDipole;
      shYOUT2->DefineSection(0, z,             rOuSaa1String,       252.);
//    Start of SAA1-SAA2
      z = ziSaa1Saa2;
      shYOUT2->DefineSection(1, z,             rOuSaa1String,       252.);
      shYOUT2->DefineSection(2, z,             rOuSaa1Saa2Steel,    252.);
//    End of SAA1-SAA2
      z = ziSaa2;
      shYOUT2->DefineSection(3, z,             rOuSaa1Saa2Steel,    252.);
//    SAA2
      shYOUT2->DefineSection( 4, z,            rInSaa2StEnv1 + dSt, 252.);
      z = ziSaa2 + zSaa2PbRing;
      shYOUT2->DefineSection( 5, z,            rInSaa2StEnv1 + dSt, 252.);
//    Pb Cone
      shYOUT2->DefineSection( 6, z,            rOuSaa2PbRingF,      252.);
      rmin = rOuSaa2PbRingF + (1380. - z) * TMath::Tan(1.6 * kDegRad);
      shYOUT2->DefineSection( 7, 1380., rmin, 252.);
      shYOUT2->DefineSection( 8, 1380., rmin, 304.);
      z = ziSaa2 + zSaa2PbRing + dzSaa2PbRing;
      shYOUT2->DefineSection( 9, z,            rOuSaa2PbRingR,      304.);
//    Straight Sections
      shYOUT2->DefineSection(10, z,            rInSaa2StEnv1 + dSt, 460.);
      z = ziSaa2 + dzSaa2StEnv1;
      shYOUT2->DefineSection(11, z,            rInSaa2StEnv1 + dSt, 460.);
      shYOUT2->DefineSection(12, z,            rInSaa2StEnv2 + dSt, 460.);
      z += dzSaa2StEnv2;
      shYOUT2->DefineSection(13, z,            rInSaa2StEnv2 + dSt, 460.);
      
      InvertPcon(shYOUT2);
      TGeoVolume* voYOUT2 = new TGeoVolume("YOUT2", shYOUT2, kMedAirMu);
      top->AddNode(voYOUT2, 1, vec0);
}

void AliSHILv3::Init()
{
  //
  // Initialise the muon shield after it has been built
  //
  Int_t i;
  //
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" SHILvTGeo_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    // Here the SHIL initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}

void AliSHILv3::InvertPcon(TGeoPcon* pcon)
{
//
//  z -> -z
// 
    Int_t nz = pcon->GetNz();
    Double_t* z    = new Double_t[nz];
    Double_t* rmin = new Double_t[nz];
    Double_t* rmax = new Double_t[nz];

    Double_t*  z0    = pcon->GetZ();
    Double_t*  rmin0 = pcon->GetRmin();
    Double_t*  rmax0 = pcon->GetRmax();
    
    for (Int_t i = 0; i < nz; i++) {
	z[i] = z0[i];
	rmin[i] = rmin0[i];
	rmax[i] = rmax0[i];	
    }
    
    for (Int_t i = 0; i < nz; i++) {
	Int_t j = nz - i - 1;
	pcon->DefineSection(i, - z[j], rmin[j], rmax[j]);
    }
    
    delete[] z;
    delete[] rmin;
    delete[] rmax;
}

TGeoPcon* AliSHILv3::MakeShapeFromTemplate(TGeoPcon* pcon, Float_t drMin, Float_t drMax)
{
    // 
    // Returns new shape based on a template changing
    // the inner radii by drMin and the outer radii by drMax.
    //
    Int_t nz = pcon->GetNz();
    TGeoPcon* cpcon = new TGeoPcon(0., 360., nz);
    for (Int_t i = 0; i < nz; i++) 
	cpcon->DefineSection(i, pcon->GetZ(i), pcon->GetRmin(i) + drMin, pcon->GetRmax(i) + drMax);
    return cpcon;
}


