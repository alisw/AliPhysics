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

/* $Id: AliAD.cxx  $ */

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                  AD (ALICE Diffractive)  Detector                     //
//                                                                       //
//  This class contains the base procedures for the AD  detector         //
//  Default geometry of 2013: 16 modules                                 //
//  All comments should be sent to :                                     //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

// --- Standard libraries ---
#include <Riostream.h>

// --- ROOT libraries ---
#include <TMath.h>
#include <TString.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoShape.h>
#include <TGeoXtru.h>
#include <TTree.h>
#include <TSystem.h>
// #include <TGeoCompositeShape.h> // included in .h
#include <TGeoGlobalMagField.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoArb8.h>
#include <TClonesArray.h>
#include <TGeoTrd2.h>
#include <TParticle.h>

#include <TH2F.h>
#include <TCanvas.h>

// --- AliRoot header files ---


#include "AliADhit.h"
#include "AliADdigit.h"
#include "AliADv1.h"
#include "AliLog.h"
#include "AliConst.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliMC.h"


ClassImp(AliADv1)
//__________________________________________________________________
AliADv1::AliADv1() : 
  AliAD(),
  fADCstruct(kTRUE),
  fADCPosition(kADCInTunnel),
  fADCLightYield(93.75),
  fADCPhotoCathodeEfficiency(0.18),
  fADALightYield(93.75),
  fADAPhotoCathodeEfficiency(0.18),
  fKeepHistory(kFALSE)
{
   // Default Constructor
    fHits = 0;
}

//_____________________________________________________________________________
AliADv1::AliADv1(const char *name, const char *title) : 
  AliAD(name,title),  
  fADCstruct(kTRUE),
  fADCPosition(kADCInTunnel),
  fADCLightYield(93.75),
  fADCPhotoCathodeEfficiency(0.18),
  fADALightYield(93.75),
  fADAPhotoCathodeEfficiency(0.18),
  fKeepHistory(kFALSE)
{
   // Standard constructor for AD Detector
  
   AliModule* pipe = gAlice->GetModule("PIPE");
   if( (!pipe) ) {
      Error("Constructor","AD needs PIPE!!!\n");
      exit(1);
   } 
   fHits = new TClonesArray("AliADhit",400);
   gAlice->GetMCApp()->AddHitList(fHits);
}

//_____________________________________________________________________________
AliADv1::~AliADv1()
{
	// default destructor
}
//_____________________________________________________________________________
void AliADv1::Init()
{
  // Initialise L3 magnet after it has been built
  Int_t i;
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" ADv1_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}

//_____________________________________________________________________________
void AliADv1::CreateGeometry()
{
  //
  // Create the geometry for the AD arrays
  //
  CreateAD();
  
}
//_____________________________________________________________________________
void AliADv1::ReadADCFromEnv() {
  if (! gSystem->Getenv("CONFIG_ADC_POS"))
    return;

  TString str(gSystem->Getenv("CONFIG_ADC_POS"));

  if (0==str.CompareTo ("kADCInTunnel")) {
    fADCPosition = kADCInTunnel;
    printf("fADCPosition set to kADCInTunnel\n");
    //
  } else if (0==str.CompareTo ("kADCInCavern")) {
    fADCPosition = kADCInCavern;
    printf("fADCPosition set to kADCInCavern\n");
    //
  } else if (0==str.CompareTo ("kADCInBoth")) {
    fADCPosition = kADCInBoth;
    printf("fADCPosition set to kADCInBoth\n");
  }
}
//_____________________________________________________________________________
TGeoCompositeShape * AliADv1::MakeShapeADCpadH(const Double_t W, const Double_t H, const Double_t dz) {
  /////////////////////////////////////////////////////////////////////////////
  ///                ADC pad in the cavern (H shapped hole)                 ///
  /////////////////////////////////////////////////////////////////////////////
  // const Double_t W = kADCCellSide; // Width  of Scintillator pad
  // const Double_t H = kADCCellSide; // Height of Scintillator pad
  // Coordinates of ADC pad vertexes
  Double_t pad0_x [] = { 10., 10.,   W,  W  };
  Double_t pad0_y [] = {  0., 11., 11.,  0. };
  Double_t pad1_x [] = { 15., 15.,   W,  W  };
  Double_t pad1_y [] = { 11., 15., 15., 11. };
  Double_t pad2_x [] = {  0.,  0.,   W,  W  };
  Double_t pad2_y [] = { 15.,  H ,   H, 15. };
  TGeoArb8 * shADCpad0H = new TGeoArb8("shADCpad0H", dz/2.);
  TGeoArb8 * shADCpad1H = new TGeoArb8("shADCpad1H", dz/2.);
  TGeoArb8 * shADCpad2H = new TGeoArb8("shADCpad2H", dz/2.);
  for (Int_t i=0; i<4; i++) {
    // -dz
    shADCpad0H -> SetVertex(i,   pad0_x[i], pad0_y[i]);
    shADCpad1H -> SetVertex(i,   pad1_x[i], pad1_y[i]);
    shADCpad2H -> SetVertex(i,   pad2_x[i], pad2_y[i]);
    // +dz
    shADCpad0H -> SetVertex(i+4, pad0_x[i], pad0_y[i]);
    shADCpad1H -> SetVertex(i+4, pad1_x[i], pad1_y[i]);
    shADCpad2H -> SetVertex(i+4, pad2_x[i], pad2_y[i]);
  }
  return new TGeoCompositeShape("shADCpadH","shADCpad0H + shADCpad1H + shADCpad2H");
}

//_____________________________________________________________________________
void AliADv1::CreateAD()
{
  printf("===> AliADv1::CreateAD(): ver=[Feb 3st, 2015]; contact=[ecalvovi@cern.ch]\n");
  //
  // Define Rotations used
  //
  TGeoRotation * Rx180, * Rz180, * Ry180, * Rx90, * Rx90m, * Ry90m;
  Rx90m = new TGeoRotation("Rx90m",   0., -90.,   0.) ;
  Rx90  = new TGeoRotation("Rx90" ,   0.,  90.,   0.) ;
  Rx180 = new TGeoRotation("Rx180",   0., 180.,   0.) ;  //   4    |   1
  Rz180 = new TGeoRotation("Rz180", 180.,   0.,   0.) ;  // --------------->  x  
  Ry180 = new TGeoRotation("Ry180", 180., 180.,   0.) ;  //   3    |   2
  Ry90m = new TGeoRotation("Ry90m",  90., -90., -90.) ;
  // Get Mediums needed.
  TGeoMedium * kMedAlu       = gGeoManager->GetMedium("AD_Alum");   // Aluminium 
  TGeoMedium * kMedSteelSh   = gGeoManager->GetMedium("AD_ST_C0");  // Stainless Steel 
  TGeoMedium * kMedVacuum    = gGeoManager->GetMedium("AD_VA_C0");  // Stainless Steel 
  //
  // Private comunication by Arturo Tauro (2014, Apr 23)
  // According to last survey measurement done this morning, 
  // the C-side wall is at Z = - 18959mm.
  //
  const Double_t kZwall        = 1895.9 ;  // Aluminium plate z position 
  const Double_t kZendAbs      = 1880.75;  // End of CC block absorber
  const Double_t kZbegVMAOI    = 1919.2 ;  // Begining of Warm Module
  const Double_t kZbegValve    = 1910.7 ;  // Begining of Valve
  const Double_t kZbegFrontBar = 1949.1 ;  // Begining of Front Bar
  const Double_t kZbegCoil     = 1959.4 ;  // Begining of compensator coil

  // Define Ion Pump ??
  // Drawing LHCVBU__0052
  // Vacuum - Bellows - U type
  // BODY 1 PORTS
  //
  (new TGeoCombiTrans("ctPumpVB2", 0., -13./2., 6.8-11.5/2., Rx90))->RegisterYourself();
  new TGeoTube("shIonPumpVB1o",  0.0, 10.3 /2., 11.5/2.   );
  new TGeoTube("shIonPumpVB2o",  0.0,  7.0 /2., 13.0/2.   );
  new TGeoTube("shIonPumpVB1i",  0.0, 10.0 /2., 11.5/2.+2.);
  new TGeoTube("shIonPumpVB2i",  0.0,  6.7 /2., 13.0/2.+2.);
  new TGeoCompositeShape("shIonPumpVBo", "shIonPumpVB1o+shIonPumpVB2o:ctPumpVB2");
  //
  // Continue definition of LHCVBU__0052
  //
  new TGeoCompositeShape("shIonPumpVBi", "shIonPumpVB1i+shIonPumpVB2i:ctPumpVB2");
  TGeoShape * sh3 = new TGeoCompositeShape("shIonPumpVB",  "shIonPumpVBo-shIonPumpVBi");
  TGeoVolume * voIonPumpVB = new TGeoVolume("voIonPumpVB", sh3, kMedAlu);
  // Variables 
  Double_t alpha, beta, tga2, tga, sa, ca, ctgb, d, Ro, Ri, phi1, dphi, H, L, z;
  // Drawing: LHCVSR__0054
  // Vacuum Screen - RF
  // transition flange
  alpha = 15. * TMath::DegToRad();
  beta  = 15. * TMath::DegToRad();
  tga2  = TMath::Tan(alpha*0.5);
  tga   = TMath::Tan(alpha);
  ctgb  = 1./TMath::Tan(beta );
  ca    = TMath::Cos(alpha    );
  d     = 0.3; // thickness of transition flange
  Double_t h = d/ca; // vertical distance between parallel surfaces tilted alpha degrees
  TGeoPcon * shVSRflange = new TGeoPcon("shVSRflange", 0.0, 360.0, 7);
  Ri = 9.71/2.; Ro = 11.16/2.;
  shVSRflange->DefineSection(0,         0.  ,   Ri, Ro);
  Ri = Ri - 0.25 * tga;
  shVSRflange->DefineSection(1,         0.25,   Ri, Ro);
  Ro = Ri + d*tga;
  shVSRflange->DefineSection(2,         0.25,   Ri, Ro);
  Ri = 6.3/2.; Ro = Ri + h;
  z = (9.71/2. - Ri) / tga;
  shVSRflange->DefineSection(3,            z,   Ri, Ro);
  // 
  Double_t   Dtga = 6.6*tga - 0.5*(9.71-6.3)    ;
  Double_t x = (h - 0.11 - Dtga) / (ctgb - tga) ;
  Double_t y = x * ctgb;
  z  = 6.6 - x;
  Ro = Ri + 0.11 + y;
  shVSRflange->DefineSection(4,            z,   Ri, Ro);
  z  = 6.6; 
  Ro = Ri + 0.11;
  shVSRflange->DefineSection(5,            z,   Ri, Ro);
  z  = 7.1; 
  shVSRflange->DefineSection(6,            z,   Ri, Ro);
  TGeoVolume * voVSRflange = new TGeoVolume("voVSRflange", shVSRflange, kMedAlu);
  //
  // Drawing: LHCVSR__0053
  // Vacuum Screen - RF
  // transition tube
  alpha = 10. * TMath::DegToRad();
  tga2  = TMath::Tan(alpha*0.5);
  tga   = TMath::Tan(alpha);
  ca    = TMath::Cos(alpha    );
  d     = 0.3; // thickness of Vacuum Screen RF
  Ro    = 6.7/2.; //
  Ri    = 0.;
  phi1  = 90. - 15.;
  dphi  = 30.;
  TGeoPcon * shVSR0 = new TGeoPcon("shVSR0", phi1, dphi, 6);
  shVSR0->DefineSection(0,         0.  ,   Ro-0.09, Ro);
  shVSR0->DefineSection(1,         0.45,   Ro-0.09, Ro);
  Ri=Ro-d;
  shVSR0->DefineSection(2,         0.45,      Ri, Ro);
  shVSR0->DefineSection(3, 13.37-d*tga2,      Ri, Ro);
  Ro  = Ri + d/ca;
  shVSR0->DefineSection(4,        13.37,      Ri, Ro);
  Ri += 0.63*tga;
  Ro  = Ri + d/ca;
  shVSR0->DefineSection(5,        14.0 ,      Ri, Ro);
  // printf("  Ro: %8.2f\n", Ro);
  // Make holes 
  new TGeoBBox("shHoleBody"    , 0.15, 0.60, 0.3);
  new TGeoTube("shHoleEnd", 0.  , 0.15, 0.3);
  (new TGeoTranslation("trHoleBody", 0., -0.6, 0.))->RegisterYourself();
  (new TGeoTranslation("trHoleEnd" , 0., -1.2, 0.))->RegisterYourself();
  new TGeoCompositeShape("shHole","shHoleEnd + shHoleEnd:trHoleEnd + shHoleBody:trHoleBody");
  // Single hole made. Now define some combitrans to position holes
  z = 1.3; Ro = (6.7 - d)*0.5;
  (new TGeoCombiTrans("ctHole1", 0., Ro , z, Rx90m))->RegisterYourself();
  z = 1.3 + 1*2.5;
  (new TGeoCombiTrans("ctHole2", 0., Ro , z, Rx90m))->RegisterYourself();
  z = 1.3 + 2*2.5;
  (new TGeoCombiTrans("ctHole3", 0., Ro , z, Rx90m))->RegisterYourself();
  z = 1.3 + 3*2.5;
  (new TGeoCombiTrans("ctHole4", 0., Ro , z, Rx90m))->RegisterYourself();
  // Now make a sector of RF transition tube
  new TGeoCompositeShape("shVSRsec",
   "shVSR0 - (shHole:ctHole1 + shHole:ctHole2 + shHole:ctHole3 + shHole:ctHole4)");
  // Now define rotations for each sector
  TString strSh = "shVSRsec ";
  for (Int_t i=1; i<=11; i++) {
   (new TGeoRotation(Form("rSec%d",i), 30. * i, 0. , 0.))->RegisterYourself();
   strSh+=Form("+ shVSRsec:rSec%d",i);
  }
  // printf("%s\n", strSh.Data());
  TGeoCompositeShape * shVSR = new TGeoCompositeShape("shVSR", strSh.Data());
  // Now assembly the sector to form VSR RF transition tube !
  TGeoVolume * voVSR = new TGeoVolume("voVSR", shVSR, kMedAlu);
  // 
  // Drawing: LHCVSR__0057
  // RF CONTACT
  // 
  Ro = 0.5 * 6.3;
  d  = 0.03;
  Ri = Ro - d;
  // alpha = TMath::ArcSin((7.35-1.75)/Ri);  <-- No!
  H = 0.5 * 6.9 - Ri;
  L = 28. - 7.1 + 0.45 -14. - 1.75;//7.35 - 1.75;
  // Double_t Delta = TMath::Sqrt( L*L + 4.*(H-d)*H );
  Double_t R = TMath::Sqrt((H-d)*(H-d) + L*L);
  alpha = TMath::ASin(d/R) + TMath::ASin((H-d)/R);
  //printf("alpha: %8.2f \n", alpha * TMath::RadToDeg());
  sa = TMath::Sin(alpha);
  ca = TMath::Cos(alpha);
  x = d*sa;
  y = d*ca;
  Double_t R0 =  1.75;
  Double_t R1 = 10.48;
  Double_t R2 =  0.81 + 0.28;
  phi1 = 0.; dphi = 360.; 

  TGeoPcon * shVSRcontact = new TGeoPcon("shVSRcont", phi1, dphi, 6);
  z = 0.;
  shVSRcontact->DefineSection(0, -z, Ri, Ro);
  z = R0;
  shVSRcontact->DefineSection(1, -z, Ri, Ro);
  z  += x;
  Ri += y;
  Ro  = Ri + d;
  shVSRcontact->DefineSection(2, -z, Ri, Ro);
  z  += (R1 - R0) * ca;
  Ri += (R1 - R0) * sa;
  Ro  = Ri + d;
  shVSRcontact->DefineSection(3, -z, Ri, Ro);
  // Last sections (R2)
  Double_t ab = alpha + 21. * TMath::DegToRad();
  Double_t sab = TMath::Sin(ab);
  Double_t cab = TMath::Cos(ab);
  x   = d * sab;
  y   = d/ca - d*cab;
  z  += x;
  Ri += y;
  Ro  = Ri + d/cab;
  shVSRcontact->DefineSection(4, -z, Ri, Ro);
  z  += R2 * cab;
  Ri += R2 * sab;
  Ro  = Ri + d/cab;
  shVSRcontact->DefineSection(5, -z, Ri, Ro);
  TGeoVolume * voVSRcontD = new TGeoVolume("voVSRcontD", shVSRcontact, kMedAlu);
  // Drawing: LHCVSR__0017
  // Vacuum Screen - RF
  // RF Contact flange
  phi1 = 0.;
  dphi = 360.;

  TGeoPcon * shVSRcontFlange = new TGeoPcon("shVSRcontFlange", phi1, dphi, 11);
  Ri = 6.30 /2.; Ro = 11.16 /2.; z = 0;
  shVSRcontFlange->DefineSection( 0, -z, Ri, Ro);
  Ri = 6.30 /2.; Ro = 11.16 /2.; z = 0.1;
  shVSRcontFlange->DefineSection( 1, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro = 11.16 /2.; z = 0.1;
  shVSRcontFlange->DefineSection( 2, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro = 11.16 /2.; z = 0.25;
  shVSRcontFlange->DefineSection( 3, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  9.30 /2.; z = 0.25;
  shVSRcontFlange->DefineSection( 4, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  9.30 /2.; z = 0.30;
  shVSRcontFlange->DefineSection( 5, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.90 /2.; z = 0.30;
  shVSRcontFlange->DefineSection( 6, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.90 /2.; z = 1.10;
  shVSRcontFlange->DefineSection( 7, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.85 /2.; z = 1.15;
  shVSRcontFlange->DefineSection( 8, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.50 /2.; z = 1.15;
  shVSRcontFlange->DefineSection( 9, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.50 /2.; z = 1.30;
  shVSRcontFlange->DefineSection(10, -z, Ri, Ro);
  TGeoVolume * voVSRcontF = new TGeoVolume("voVSRcontF", shVSRcontFlange, kMedAlu);
  // Drawing: LHCVBU__0002
  // Bellows + End Parts
  // Vacuum - Bellows - U type

  // First make end part
  phi1 = 0. ; dphi = 360. ;
  TGeoPcon * shVBUend = new TGeoPcon("shVBUend", phi1, dphi, 6);
  Ri = 5.176; Ro = 5.4; z = 0;
  shVBUend->DefineSection( 0, z, Ri, Ro);
  Double_t dz  = 0.03;
  Ri -= dz * TMath::Tan(45. * TMath::DegToRad());
  z  += dz;
  shVBUend->DefineSection( 1, z, Ri, Ro);
  dz  = 0.04;
  Ro -= dz * TMath::Tan(15. * TMath::DegToRad());
  Ri -= dz * TMath::Tan(45. * TMath::DegToRad());
  z  += dz;
  shVBUend->DefineSection( 2, z, Ri, Ro);
  Ro  = 5.250 ;
  Ri  = 5.073;
  z   = 0.103;
  shVBUend->DefineSection( 3, z, Ri, Ro);
  Ro  = 5.25;
  Ri  = 5.02;
  z   = 0.18;
  shVBUend->DefineSection( 4, z, Ri, Ro);
  Ro  = 5.15;
  Ri  = 5.00;
  z   = 0.28;
  shVBUend->DefineSection( 5, z, Ri, Ro);

  TGeoPcon * shVBUtube26mm = new TGeoPcon("shVBUtube26mm", 0., 360., 2);
  shVBUtube26mm->DefineSection( 0, 0.28, 5.0, 5.15);
  shVBUtube26mm->DefineSection( 1, 2.60, 5.0, 5.15);
  TGeoPcon * shVBUtube9mm  = new TGeoPcon("shVBUtube9mm" , 0., 360., 2);
  shVBUtube9mm ->DefineSection( 0, 0.00, 5.0, 5.15);
  shVBUtube9mm ->DefineSection( 1, 0.62, 5.0, 5.15);
  // central part of bellow (TODO: Add plies)
  TGeoPcon * shVBUcent  = new TGeoPcon("shVBUcent" , 0., 360., 6);
  const Int_t nsec = 6;
  Double_t az [nsec] = {0.9, 0.915, 0.915, 11.885, 11.885, 11.9};
  Double_t aRi[nsec] = {5.0, 5.0  , 5.685,  5.685,  5.0  ,  5.0};
  for (Int_t i=0; i<nsec; i++) {
   z=az[i]; Ri = aRi[i]; Ro = 5.7;
   //printf("  i: %2d  z: %8.2f  Ri: %8.2f  Ro: %8.2f\n", i, z, Ri, Ro );
   shVBUcent ->DefineSection( i, z, Ri, Ro);
  }

  ( new TGeoCombiTrans("ctEnd9mm", 0., 0., 0.9, Ry180)) -> RegisterYourself();

  TGeoCompositeShape * shVBU9mm  = new TGeoCompositeShape("shVBU9mm" , "shVBUend:ctEnd9mm + shVBUtube9mm");
  TGeoCompositeShape * shVBU26mm = new TGeoCompositeShape("shVBU26mm", "shVBUend + shVBUtube26mm");

  const Int_t nsec2 = 10;
  Double_t az2 [nsec2] = {0.  , 0.43, 0.43, 0.58, 0.58 , 0.73 , 0.73, 1.05, 1.05, 1.31} ;
  Double_t aRi2[nsec2] = {5.15, 5.15, 5.03, 5.03, 5.455, 5.455, 5.03, 5.03, 5.59, 5.59} ;
  TGeoPcon * shVBUrotFlg  = new TGeoPcon("shVBUrotFlg" , 0., 360., nsec2);
  for (Int_t i=0; i<nsec2; i++) {
   z=az2[i]; Ri = aRi2[i]; Ro = 6.02;
   //printf("  i: %2d  z: %8.2f  Ri: %8.2f  Ro: %8.2f\n", i, z, Ri, Ro );
   shVBUrotFlg ->DefineSection( i, z, Ri, Ro);
  }
  TGeoVolume * voVBUrotFlg = new TGeoVolume("voVBUrotFlg", shVBUrotFlg, kMedAlu);

  // Flange 
  TGeoPcon * shVBUflg  = new TGeoPcon("shVBUflg" , 0., 360., 4);
  z  = 0; 
  Ri = 6.02; Ro = 7.6; 
  shVBUflg->DefineSection(0, z, Ri, Ro);
  z  = 1.31; 
  shVBUflg->DefineSection(1, z, Ri, Ro);
  z  = 1.31; 
  Ri = 5.15; 
  shVBUflg->DefineSection(2, z, Ri, Ro);
  z  = 2.0;
  shVBUflg->DefineSection(3, z, Ri, Ro);
  TGeoVolume * voVBUflg = new TGeoVolume("voVBUflg", shVBUflg, kMedAlu);

   
  // Add the metal plate at the end of Absorber
  // Plate:  80 cm x 80 cm x 1.95 cm (Thickness is aproximatted) (ernesto.calvo@pucp.edu.pe)
  // The End of the concrete abosorber is at kZendAbs = 1880.75

  new TGeoBBox("shBasePlate", 80./2., 80./2.,  1.95/2.);
  new TGeoTube("shHolePlate",  0.   , 12.3  ,  1.95   );
  TGeoVolume* voSaa3EndPlate  =  new TGeoVolume("voYSaa3EndPlate",
      new TGeoCompositeShape("shYSaa3EndPlate","shBasePlate-shHolePlate"),
      kMedSteelSh);
  //
  // Add Rods
  //
  // dimensions of rods
  const Double_t dzRodL = 27.0;
  const Double_t dzRodA =  4.3;
  const Double_t dzRodB =  1.3;
  new TGeoTube("shLargeRod", 0.,   1.6/2.,  dzRodL/2.);
  new TGeoTube("shRodA",     0.,   3.0/2.,  dzRodA/2.);
  new TGeoTube("shRodB",     0.,   2.3/2.,  dzRodB/2.);
  //
  ( new TGeoTranslation("trRod1", 0., 0., -dzRodL/2. + dzRodA/2.)          )->RegisterYourself();
  ( new TGeoTranslation("trRod2", 0., 0., -dzRodL/2. + dzRodA + dzRodB/2.) )->RegisterYourself();
  ( new TGeoTranslation("trRod3", 0., 0.,  dzRodL/2. - dzRodB/2.)          )->RegisterYourself();
  TGeoVolume * voSaa3Rod = new TGeoVolume("YSAA3Rod",
      new TGeoCompositeShape("shLargeRod+shRodA:trRod1 + shRodB:trRod2 + shRodB:trRod3"),
      kMedSteelSh);
  //
  // Define Valve support (VS)  (ernesto.calvo@pucp.edu.pe) 
  //
  const Double_t dyVS =  5.5; 
  const Double_t dxVS = 30.0; 
  const Double_t dzVS =  1.0;
  TGeoVolume * voVS = new TGeoVolume("voVS",
      new TGeoBBox("shVS", dxVS/2., dyVS/2., dzVS/2.),
      kMedSteelSh);

  // Add Valve (Valve is divided in parts VA,VB,VC and VD)  (ernesto.calvo@pucp.edu.pe) 
  TGeoVolumeAssembly * voValve = new TGeoVolumeAssembly("voValve");
  //
  // Define volume VA  (ernesto.calvo@pucp.edu.pe) 
  const Double_t dxVA  = 20.3; 
  const Double_t dyVA  = 48.0; 
  const Double_t dzVA  =  6.0; // Width 
  const Double_t dz2VA =  8.5; // Full width including protuding boxes
  // Valve position  (ernesto.calvo@pucp.edu.pe) 
  const Double_t zPosValve =  kZbegValve + dz2VA/2.;
  //
  new TGeoBBox("shVAbox",       dxVA/2., dyVA/2.,     dzVA/2.);
  new TGeoBBox("shVAHbox",  -1.+dxVA/2.,   3./2.,    dz2VA/2.);
  new TGeoTube("shVAC",              0.,     7.9,    dz2VA/2.);
  new TGeoTube("shVACh",             0.,     5.3,    dz2VA   );
  // translation for shVAHbox (ernesto.calvo@pucp.edu.pe) 
  ( new TGeoTranslation("trVAH1", 0.,  12.75, 0.) )->RegisterYourself();
  ( new TGeoTranslation("trVAH2", 0., -12.75, 0.) )->RegisterYourself();

  TGeoVolume * voValveVA = new TGeoVolume("voValveVA",
      new TGeoCompositeShape("(shVAbox + shVAHbox:trVAH1 + shVAHbox:trVAH2 + shVAC)-shVACh"),
      kMedSteelSh);
  voValve->AddNode(voValveVA, 1, 0);
  // Define Vacuum Hole of Valve
  TGeoTube   * shVACvach   = new TGeoTube("shVACvach", 0., 5.3, dz2VA/2.);
  TGeoVolume * voValveVAvh = new TGeoVolume("voValveVAvacuum", shVACvach, kMedVacuum);
  voValve->AddNode(voValveVAvh,1,0);
  // Also add valve Support (ernesto.calvo@pucp.edu.pe) 
  voValve->AddNode(voVS, 1, new TGeoTranslation(0.,  12.75, -dz2VA/2.-dzVS/2.));
  voValve->AddNode(voVS, 2, new TGeoTranslation(0., -12.75, -dz2VA/2.-dzVS/2.));

  // Define Volume VB (ernesto.calvo@pucp.edu.pe) 
  const Double_t dxVB = 23.5; 
  const Double_t dyVB =  5.0; 
  const Double_t dzVB =  9.4; 
  TGeoVolume * voValveVB = new TGeoVolume("voValveVB", 
      new TGeoBBox("shVBbox", dxVB/2., dyVB/2., dzVB/2.),
      kMedSteelSh);
  voValve->AddNode(voValveVB, 1, new TGeoTranslation(  0., dyVA/2. +dyVB/2. , 0));
  // Define Volume VC (ernesto.calvo@pucp.edu.pe) 
  const Double_t R1VC  =  4.5 /2.;
  const Double_t R2VC  =  8.1 /2.;
  const Double_t dy1VC  = 10.0;
  const Double_t dy2VC =  0.75;
  new TGeoTube("shVC1",      0.,   R1VC, dy1VC/2.);
  new TGeoTube("shVC2",      0.,   R2VC, dy2VC/2.);
  ( new TGeoTranslation("trVC21", 0., 0.,  dy1VC/2. - dy2VC/2.) )->RegisterYourself();
  ( new TGeoTranslation("trVC22", 0., 0., -dy1VC/2. + dy2VC/2.) )->RegisterYourself();
  TGeoVolume * voValveVC = new TGeoVolume("voValveVC",
      new TGeoCompositeShape("shVC1  + shVC2:trVC21 + shVC2:trVC22"),
      kMedSteelSh);
  voValve->AddNode(voValveVC, 1, new TGeoCombiTrans( 
        0., dyVA/2. + dyVB + dy1VC/2. , 0, Rx90) );
  // Define volume VD (ernesto.calvo@pucp.edu.pe) 
  const Double_t dxVD = 15.9;
  const Double_t dyVD = 23.0;
  const Double_t dzVD = 14.0;
  TGeoVolume * voValveVD = new TGeoVolume("voValveVD",
      new TGeoBBox("shVD", dxVD/2., dyVD/2., dzVD/2.),
      kMedSteelSh);
  voValve->AddNode(voValveVD, 1, 
      new TGeoTranslation( 1.25, dyVA/2. + dyVB + dy1VC + dyVD/2. , 0) );

  //
  // Define volume Front Bar (ernesto.calvo@pucp.edu.pe) 
  //
  const Double_t dxF  = 67.4; 
  const Double_t dyF  =  4.0; 
  const Double_t dzF  =  2.0; 
  const Double_t R1FC =  8.1;
  const Double_t R2FC = 11.5;
  const Double_t dxFA  = (dxF-R1FC)/2.; 
 
  new TGeoBBox("shFA", dxFA/2., dyF/2., dzF/2.);
  new TGeoTube("shFC",      0.,   R2FC, dzF/2.);
  new TGeoTube("shFCH",     0.,   R1FC, 2.*dzF);
  ( new TGeoTranslation("trFA1",  R1FC +dxFA/2., 0., 0.) )->RegisterYourself();
  ( new TGeoTranslation("trFA2", -R1FC -dxFA/2., 0., 0.) )->RegisterYourself();
  TGeoVolume * voFrontBar = new TGeoVolume("voFrontBar",
      new TGeoCompositeShape("shFA:trFA1 + shFA:trFA2 + (shFC - shFCH)"),
      kMedSteelSh);
  // Make Lateral Bars
  const Double_t kdzLatBar = 22.9;
  TGeoVolume * voLatBar = new TGeoVolume("voLatBar", 
    new TGeoTube("shLatBar",0., dyF/2., kdzLatBar/2.), kMedSteelSh);
  //
  // Define Compensator Magnet coils
  //
  dz = 12.5;
  Ro = 15.0;
  Ri = 1.9;
  new TGeoTubeSeg("shCoilRo", 0., Ro, dz/2., 90., 180.);
  new TGeoTubeSeg("shCoilRi", 0., Ri, dz   , 90., 185.);
  (new TGeoTranslation("trCoilRo", Ro, 10.4, 0.)) -> RegisterYourself();
  (new TGeoTranslation("trCoilRi", Ro,  9.6, 0.)) -> RegisterYourself();
  new TGeoBBox("shBoxCoilRo", 15.0/2., 10.4/2., dz/2.);
  new TGeoBBox("shBoxCoilRi",  1.9/2.,  9.6   , dz   );
  (new TGeoTranslation("trBoxCoilRo",       15.0/2., 10.4/2., 0.)) -> RegisterYourself();
  (new TGeoTranslation("trBoxCoilRi", 13.1 + 1.9/2.,  0.    , 0.)) -> RegisterYourself();
  new TGeoBBox("shBoxCoil0", 10./2., 30., dz);
  (new TGeoTranslation("trBoxCoil0",  14.6 + 10./2.,  0.    , 0.)) -> RegisterYourself();
  strSh  = "";
  strSh += "(shCoilRo:trCoilRo + shBoxCoilRo:trBoxCoilRo) - ";
  strSh += "(shCoilRi:trCoilRi + shBoxCoilRi:trBoxCoilRi +"  ;
  strSh += " shBoxCoil0:trBoxCoil0 )"  ;
  TGeoCompositeShape * shCoil = new TGeoCompositeShape("shCoil0", strSh);
  TGeoVolume * voCoil = new TGeoVolume("voCoil", shCoil, kMedSteelSh);

  // 
  // ALUMINIUM PLATES 
  //

  // Shape for aluminium Plate separating cavern and LHC tunnel
  const Double_t dAlWallThick = 0.5; // thickness of aluminium plates (cm)
  //
  // RB24/26 Tunnel Floor 
  R   = 220.;
  // h   = 140.;
  // phi = TMath::ACos(h / r);
  // xl  = r * TMath::Sin(phi);
  // dr  = 1600.;
  // dh  = dr * TMath::Cos(phi);
  // dl  = dr * TMath::Sin(phi);

  new TGeoTube("shWallBase",    0.,    R, dAlWallThick*0.5);  // base shape for shWallBigPlate
  new TGeoBBox("shWallCutBot",270., 110., dAlWallThick    );  // to be substracted from base
  // Translation for cutting circular and square hole in the plates
  (new TGeoTranslation("trAntiBeamAxis",   -70.,               40., 0.)) -> RegisterYourself();
  (new TGeoTranslation("trHUWAT3",         -70., -110. - 140. +40., 0.)) -> RegisterYourself();

  //
  //  Wall Big Aluminium Plate with Squared Hole 
  //
  const Double_t dSqrHoleSide = 33.0; // Side
  new TGeoBBox("shWallSqrHole", dSqrHoleSide*0.5, dSqrHoleSide*0.5, dAlWallThick);
  strSh  = ""; 
  strSh += "shWallBase:trAntiBeamAxis - ";
  strSh += " ( shWallCutBot:trHUWAT3" ;
  strSh += " + shWallSqrHole )";
  TGeoVolume* voWallBigPlate = new TGeoVolume("voWallBigPlate", 
    new TGeoCompositeShape("shWallBigPlate", strSh), kMedAlu );
  //
  // Wall Squared Aluminium Plate 
  //
  const Double_t dCircHoleRad = 9.5; // Radius
  new TGeoTube("shCircHole", 0., dCircHoleRad, dAlWallThick);
  // Make holes for bars
  new TGeoTube("shRodHole", 0.,   1.7/2.,  2*dAlWallThick);
  ( new TGeoTranslation("trWallRod1",  12.5, -12.75, 0.)) -> RegisterYourself();
  ( new TGeoTranslation("trWallRod2",  12.5,  12.75, 0.)) -> RegisterYourself();
  ( new TGeoTranslation("trWallRod3", -12.5, -12.75, 0.)) -> RegisterYourself();
  ( new TGeoTranslation("trWallRod4", -12.5,  12.75, 0.)) -> RegisterYourself();
  new TGeoBBox("shWallSqrPlateBase", dSqrHoleSide*0.5 + 5.0, dSqrHoleSide*0.5 + 5.0, dAlWallThick/2.);
  strSh  = " ((((";
  strSh += " ( shWallSqrPlateBase - shCircHole )"; 
  strSh += " - shRodHole:trWallRod1)"  ;
  strSh += " - shRodHole:trWallRod2)"  ;
  strSh += " - shRodHole:trWallRod3)"  ;
  strSh += " - shRodHole:trWallRod4)"  ;
  TGeoVolume* voWallSqrPlate = new TGeoVolume("HUWAT_AlWall02", 
    new TGeoCompositeShape("shWallSqrPlate", strSh ), kMedAlu);
  // ==========================================================================
  //
  // Define Mother Vacuum volume of VMAOI  (need shIonPumpVBo)
  //
  // ==========================================================================
  const Double_t kdzMoFlange   =  2.0;
  const Double_t kdzMoBellow   = 15.6;
  const Double_t kdzTTube = 11.5; // Bellow starts here
  const Double_t kziTTube  =  1.0; // Ion Pum Tube starts here
  new TGeoTube( "shMoFlange", 0., 15.2/2., kdzMoFlange/2.0);
  new TGeoTube( "shMoBellow", 0., 11.4/2., kdzMoBellow/2.0);
  (new TGeoTranslation("trMoFlange1", 0., 0.,        0.5*kdzMoFlange     )) -> RegisterYourself();
  (new TGeoTranslation("trMoFlange2", 0., 0., 28.0 - 0.5*kdzMoFlange     )) -> RegisterYourself();
  (new TGeoTranslation("trMoBellow" , 0., 0., 28.0 - 0.5*kdzMoBellow     )) -> RegisterYourself();
  (new TGeoTranslation("trMoTTube"  , 0., 0., kziTTube +  0.5*kdzTTube   )) -> RegisterYourself();
  
  TGeoCompositeShape * shMoVMAOI =  
      new TGeoCompositeShape("shMoVMAOI", 
      "shMoFlange:trMoFlange1 + shMoBellow:trMoBellow + shMoFlange:trMoFlange2 + shIonPumpVBo:trMoTTube");
  TGeoVolume * voMoVMAOI = new TGeoVolume("voMoVMAOI", shMoVMAOI, kMedVacuum);
  voMoVMAOI->AddNode(voVSR      ,1, new TGeoTranslation(0., 0., 7.1 - 0.45));
  voMoVMAOI->AddNode(voIonPumpVB,1, new TGeoTranslation(0., 0., 1 + 11.5/2.));
  voMoVMAOI->AddNode(voVSRflange,1, new TGeoTranslation(0.,0.,0.));
  voMoVMAOI->AddNode(voVSRcontD ,1, new TGeoTranslation(0.,0.,28.));
  voMoVMAOI->AddNode(voVSRcontF ,1, new TGeoTranslation(0.,0.,28.));
  z = 1.0 + 11.5;
  voMoVMAOI->AddNode( new TGeoVolume("voVBU9mm", shVBU9mm, kMedAlu),
               1, new TGeoTranslation(0.,0., z) );
  voMoVMAOI->AddNode( new TGeoVolume("voVBU26mm", shVBU26mm, kMedAlu),
               1, new TGeoTranslation(0.,0., z + 11.9) );
  voMoVMAOI->AddNode( new TGeoVolume("voVBUcent", shVBUcent, kMedAlu),
               1, new TGeoTranslation(0.,0., z) );
  voMoVMAOI->AddNode( voVBUrotFlg, 1, 
               new TGeoCombiTrans(0.,0.,1.31, Ry180) );
  voMoVMAOI->AddNode( voVBUrotFlg, 2, 
               new TGeoTranslation(0.,0.,28. - 1.31) );
  voMoVMAOI->AddNode( voVBUflg, 1, 
               new TGeoTranslation(0.,0.,0.) );
  voMoVMAOI->AddNode( voVBUflg, 2, 
               new TGeoCombiTrans(0.,0.,28., Ry180) );
  // ==========================================================================
  //
  // AD Support structure by Pieter Ijzerman
  // ecalvovi@cern.ch
  // ==========================================================================
  Int_t nvertices=0;
  // Cover plate_______________________________________________________________
  TGeoXtru * shADcoverplate = new TGeoXtru(2);
  shADcoverplate->SetNameTitle("shADcoverplate","shADcoverplate");
  Double_t y1[] = {  0.0, 18.50, 18.50, 22.50, 22.50, 18.50, 18.50, 22.50, 22.50, 18.50, 18.50,   .00 ,  .00, 15.25, 15.25,  .00 }; 
  Double_t x1[] = {  0.0,   .00,  5.15,  5.15, 17.15, 17.15, 24.25, 24.25, 36.25, 36.25, 41.40, 41.40 ,35.70, 35.70,  5.70, 5.70 }; 
  nvertices = sizeof(x1)/sizeof(Double_t);
  shADcoverplate->DefinePolygon(nvertices,x1,y1);
  shADcoverplate->DefineSection(0, -0.1, -20.7, 0.0, 1.0); // Z position, offset and scale for first section
  shADcoverplate->DefineSection(1,  0.1, -20.7, 0.0, 1.0); // -''- secons section

  // Horizontal side___________________________________________________________
  TGeoXtru * shADhorizontalside = new TGeoXtru(2);
  shADhorizontalside->SetNameTitle("shADhorizontalside","shADhorizontalside");
  Double_t x2[] = {  0.0,  .00, 4.80, 4.80, 7.20, 7.20, 12.00, 12.00 };
  Double_t y2[] = {  0.0, 5.66, 5.66, 1.16, 1.16, 5.66,  5.66,   .00 };
  nvertices = sizeof(x2)/sizeof(Double_t);
  shADhorizontalside->DefinePolygon(nvertices,x2,y2);
  shADhorizontalside->DefineSection(0, -0.4, -6.0, 0.0, 1.0); // Z position, offset and scale for first section
  shADhorizontalside->DefineSection(1, +0.4, -6.0, 0.0, 1.0); // -''- secons section

  TGeoBBox * shADsidebox = new TGeoBBox("shADsidebox", 0.4, 18.55/2., 5.66/2.);
  TGeoVolume * voADsidebox = new TGeoVolume("voADsidebox", shADsidebox, kMedAlu);


// Define a TNode where this example resides in the TGeometry
// Draw the TGeometry
  TGeoVolume * voADhorizontalside = new TGeoVolume("voADhorizontalside", shADhorizontalside, kMedAlu);
  TGeoVolume * voADcoverplate = new TGeoVolume("voADcoverplate", shADcoverplate, kMedAlu);
  //
  TGeoVolume *voADsupport = new TGeoVolumeAssembly("voADsupport"); 
  voADsupport->AddNode(voADcoverplate,  1, new TGeoTranslation( 0., 0., -5.66/2.-0.23));
  voADsupport->AddNode(voADcoverplate,  2, new TGeoTranslation( 0., 0., +5.66/2.+0.23));
  voADsupport->AddNode(voADhorizontalside,  1, new TGeoCombiTrans( -6.0 - 7.1/2., 22.5-0.4, -5.66/2., Rx90));
  voADsupport->AddNode(voADhorizontalside,  2, new TGeoCombiTrans( +6.0 + 7.1/2., 22.5-0.4, -5.66/2., Rx90));
  voADsupport->AddNode(voADsidebox,  1, new TGeoTranslation( -20.7 +0.4, 18.55/2., 0.));
  voADsupport->AddNode(voADsidebox,  2, new TGeoTranslation( +20.7 -0.4, 18.55/2., 0.));

  // ==========================================================================
  //
  // Define ADA
  //
  // ==========================================================================

  TGeoVolume *ad = new TGeoVolumeAssembly("AD");
  
  // Get medium for ADA
  TGeoMedium * medADASci        = gGeoManager->GetMedium("AD_BC404"); // AD Scin.
  // TGeoMedium * medADALG      = gGeoManager->GetMedium("AD_PMMA");  // lightGuide
  // TGeoMedium * medADAPMGlass = gGeoManager->GetMedium("AD_Glass"); // Glass for Aluminium simulation
  // TGeoMedium * medADAPMAlum  = gGeoManager->GetMedium("AD_Alum");
  
  // Get Medium for ADC 
  TGeoMedium * medADCSci     = gGeoManager->GetMedium("AD_BC404");
  // TGeoMedium * medADCLG      = gGeoManager->GetMedium("AD_PMMA");
  // TGeoMedium * medADCPMGlass = gGeoManager->GetMedium("AD_Glass");
  // TGeoMedium * medADCPMAlum  = gGeoManager->GetMedium("AD_Alum");
  
  // ADA Scintillator Pad 
  const Double_t kADACellSideY = 21.6;
  const Double_t kADACellSideX = 18.1;
  // ADC Scintillator Pad 
  const Double_t kADCCellSideY = 21.6;
  const Double_t kADCCellSideX = 18.1;
  // WLS bar          :  0.40 cm ( 4.0 mm )
  // Wrapping         :  0.20 cm ( 2.0 mm )
  // Aluminnized Mylar:  0.01 cm ( 0.1 mm )
  // Fishing line     :  0.04 cm ( 0.4 mm )
  // total shift on X :  0.65 cm
  // total shift on Y :  0.21 cm
  const Double_t kShiftX       =  0.54;
  const Double_t kShiftY       =  0.10;
  const Double_t kADACelldz    =  2.54;
  const Double_t kADCCelldz    =  2.54;
  const Double_t kADABeamPipeR =  6.20; // Radius of beam pipe hole for ADA (Diameter  12.4 cm)
  const Double_t kADCBeamPipeR =  3.70; // Radius of beam pipe hole for ADC (Diameter   7.4 cm)
  const Int_t    kColorADA     = kGreen;
  const Int_t    kColorADC     = kGreen;
  Double_t X = kShiftX + kADACellSideX * 0.5;
  Double_t Y = kShiftY + kADACellSideY * 0.5;
  Double_t WLS_dx =  0.4;
  Double_t WLS_dz =  2.5;
  Double_t WLS_SideA_Long_dy  = 24.20; // 24.2;
  Double_t WLS_SideC_Long_dy  = 24.20; // 24.2;
  Double_t WLS_SideA_Short_dy = 18.20; // 18.41; 
  Double_t WLS_SideC_Short_dy = 20.70; // 20.91; 
  // Creating ADA WLS bars_____________________________________________________
  TGeoVolume * vADA_WLS_s = new TGeoVolume( "ADAWLSshort", 
      new TGeoBBox( "shADAWLSbarShort" , WLS_dx/2.0, WLS_SideA_Short_dy/2.0, WLS_dz/2.0),
      medADASci);      
  TGeoVolume * vADA_WLS_l = new TGeoVolume( "ADAWLSlong" , 
      new TGeoBBox( "shADAWLSbarLong"  , WLS_dx/2.0, WLS_SideA_Long_dy /2.0, WLS_dz/2.0),
      medADASci);      
  vADA_WLS_l->SetLineColor( kRed );
  vADA_WLS_s->SetLineColor( kRed );
  // Creating ADC WLS bars_____________________________________________________
  TGeoVolume * vADC_WLS_s = new TGeoVolume( "ADCWLSshort", 
      new TGeoBBox( "shADCWLSbarShort" , WLS_dx/2.0, WLS_SideC_Short_dy/2.0, WLS_dz/2.0),
      medADCSci);      
  TGeoVolume * vADC_WLS_l = new TGeoVolume( "ADCWLSlong" , 
      new TGeoBBox( "shADCWLSbarLong"  , WLS_dx/2.0, WLS_SideC_Long_dy /2.0, WLS_dz/2.0),
      medADCSci);      
  vADC_WLS_l->SetLineColor(kRed);
  vADC_WLS_s->SetLineColor(kRed);
  // Make ADA scintillator pad_________________________________________________
  new TGeoBBox( "shADAbox" , kADACellSideX/2.0, kADACellSideY/2.0, kADACelldz/2.0 );
  new TGeoTube( "shADAHole",               0. , kADABeamPipeR    , kADACelldz     );
  ( new TGeoTranslation("trADAbox", X, Y, 0.)) -> RegisterYourself();
  // 
  TGeoVolume * vADA1 = new TGeoVolume( "ADApad", 
    new TGeoCompositeShape("shADApad", "shADAbox:trADAbox-shADAHole"), medADASci );      
  vADA1->SetLineColor( kColorADA );
  
  TGeoVolume *secADA  = new TGeoVolumeAssembly( "ADAsec" ); 
  // Add PAD
  secADA->AddNode( vADA1, 1, 0); 
  secADA->AddNode( vADA_WLS_s, 1, 
      new TGeoTranslation(0.1 + WLS_dx/2.0, kADABeamPipeR + WLS_SideA_Short_dy/2.0, 0.0) ); 
  secADA->AddNode( vADA_WLS_l, 1, 
      new TGeoTranslation(kShiftX + WLS_dx/2.0 + kADACellSideX + 0.04, kShiftY + WLS_SideA_Long_dy/2.0, 0.0) ); 

  /// Assembling ADA adding 4 sectors                                       //  Sectors
  TGeoVolume *vADAarray = new TGeoVolumeAssembly( "ADA" );                  //        ^ y
  vADAarray->AddNode( secADA, 1 );                                          //        |   
  vADAarray->AddNode( secADA, 2, Ry180 );                                   //   2    |   1
  vADAarray->AddNode( secADA, 3, Rz180 );                                   // --------------->  x     
  vADAarray->AddNode( secADA, 4, Rx180 );                                   //   3    |   4
  //                                                                        //        |
  // Add ADA layer 2 and 3 to AD volume
  // const Float_t kPosAD2 = 1695.0;
  // const Float_t kPosAD3 = 1700.0;
  // ad->AddNode(vADAarray,1, new TGeoTranslation(0., 0., kPosAD2)); 
  // ad->AddNode(vADAarray,2, new TGeoTranslation(0., 0., kPosAD3));
  // const Float_t kPosADA = 1699.7;  // z-center of assembly (cm) Old
  const Float_t kPosADA = 1696.67;  // z-center of assembly (cm) New, according to Survey by F. Klumb and E.Calvo 2015 Sept 4th.
  ad->AddNode(vADAarray,    1, new TGeoTranslation(0., 0., kPosADA - kADACelldz/2. -0.23)); 
  ad->AddNode(vADAarray,    2, new TGeoTranslation(0., 0., kPosADA + kADACelldz/2. +0.23));
  ad->AddNode(voADsupport,  1, new TGeoTranslation(0., 0., kPosADA));
  ad->AddNode(voADsupport,  2, new TGeoCombiTrans (0., 0., kPosADA, Rz180));

  // ==========================================================================
  //
  // Define ADC (2014, May 4) Updated 2015, Jan 22
  //
  // ==========================================================================

  /////////////////////////////////////////////////////////////////////////////
  /// ADC in the tunnel                                                     ///
  /////////////////////////////////////////////////////////////////////////////
  new TGeoBBox( "shADCbox" , kADCCellSideX/2.0, kADCCellSideY/2.0, kADCCelldz/2.0 );
  new TGeoTube( "shADCHole",               0. , kADCBeamPipeR    , kADCCelldz     );
  X = kShiftX + kADCCellSideX * 0.5;
  Y = kShiftY + kADCCellSideY * 0.5;
  ( new TGeoTranslation("trADCbox", X, Y, 0.) ) -> RegisterYourself();
  // 
  TGeoVolume * vADCpad = new TGeoVolume( "ADCpad", 
    new TGeoCompositeShape("shADCpad", "shADCbox:trADCbox-shADCHole"), medADCSci );      
  vADCpad->SetLineColor( kColorADC );
  
  /// Creating Sector for Tunnel (Asembly:  Scintillator Pad + Light guide + PM )
  TGeoVolume *voADC  = new TGeoVolumeAssembly("ADCsec");
  // Add PAD
  voADC->AddNode( vADCpad, 1, 0);
  // Add ADC WLS Short bar
  voADC->AddNode( vADC_WLS_s, 1, 
      new TGeoTranslation( 0.1 + WLS_dx/2.0, kADCBeamPipeR + WLS_SideC_Short_dy/2.0, 0.0) ); 
  // Add ADC WLS Long  bar
  voADC->AddNode( vADC_WLS_l, 1, 
      new TGeoTranslation( 0.04 + WLS_dx/2.0 + kADCCellSideX + kShiftX, kShiftY + WLS_SideC_Long_dy/2.0, 0.0) ); 
  
  /// Assembling ADC adding the 4 sectors                 //  Sectors
  TGeoVolume *vADCarray = new TGeoVolumeAssembly("ADC");  //        ^ y
  vADCarray->AddNode( voADC, 1 );                         //        |   
  vADCarray->AddNode( voADC, 2, Ry180 );                  //   2    |   1
  vADCarray->AddNode( voADC, 3, Rz180 );                  // --------------->  x  
  vADCarray->AddNode( voADC, 4, Rx180 );                  //   3    |   4
                                                          //        |
                                                                             

  // ==========================================================================
  //
  // Add ADC to AD volume
  //
  // Note to future maintainers: 
  // In previous AliRoot versions the position z = -1900.75 corresponded 
  // to the end of the YSAA3_CC_BLOCK (concrete block shielding just before 
  // the C-Side LHC wall). Now this has been fixed to agree with reality. 
  // The YSAA3_CC_BLOCK starts at 1800.75 and ends at 1880.75 cm.
  //
  // Ernesto Calvo and Alberto Gago.
  // - ecalvovi@cern.ch
  // - agago@pucp.edu.pe
  //
  // ==========================================================================
  
  // *  ad -> AddNode(vADCarray , 1, new TGeoTranslation(0., 0., kZendADC2 + kADCCelldz/2.)); // Tunnel
  // *  ad -> AddNode(vADCarray , 2, new TGeoTranslation(0., 0., kZbegADC1 - kADCCelldz/2.)); // Tunnel
  // *  const Float_t kZbegADC1 = -kZbegFrontBar-2.
  // *  const Float_t kZendADC2 = -1959.0;            // (ecalvovi@cern.ch) 
  
  switch (fADCPosition ) {
    case kADCInTunnel:
      {
        // const Float_t kZbegADC1 = -kZbegFrontBar-2.;  // (ecalvovi@cern.ch) 
        // const Float_t kZendADC2 = -1959.0;            // (ecalvovi@cern.ch) 
        // ad -> AddNode(vADCarray , 1, new TGeoTranslation(0., 0., kZendADC2 + kADCCelldz/2.)); // Tunnel
        // ad -> AddNode(vADCarray , 2, new TGeoTranslation(0., 0., kZbegADC1 - kADCCelldz/2.)); // Tunnel
        const Float_t kPosADC = -kZbegFrontBar-2.-3.0-0.3;  // 3.0 = (5.6 + 0.2 + 0.2)/2. // (ecalvovi@cern.ch) 
        printf("CreateAD: kPosADC=%8.2f\n", kPosADC);
        ad -> AddNode(vADCarray,   1, new TGeoTranslation(0., 0., kPosADC - kADCCelldz/2. - 0.23)); // Tunnel // ADC1
        ad -> AddNode(vADCarray,   2, new TGeoTranslation(0., 0., kPosADC + kADCCelldz/2. + 0.23)); // Tunnel // ADC2
        ad -> AddNode(voADsupport, 3, new TGeoTranslation(0., 0., kPosADC));
        ad -> AddNode(voADsupport, 4, new TGeoCombiTrans (0., 0., kPosADC, Rz180));
        break;
      }
    case kADCInCavern:
      {
        printf("FATAL: vADCInCavern is now obsolete!");
        exit(1);
        // const Float_t kZbegADC1 = -1890.0;  // (ecalvovi@cern.ch) 
        // const Float_t kZbegADC2 = -1885.0;  // (ecalvovi@cern.ch) 
        // ad -> AddNode(vADCarrayH, 1, new TGeoTranslation(0., 0., kZbegADC1 - kADCCelldz/2.)); // Cavern
        // ad -> AddNode(vADCarrayH, 2, new TGeoTranslation(0., 0., kZbegADC2 - kADCCelldz/2.)); // Cavern
        break;
      }
    case kADCInBoth:
      {
        printf("FATAL: vADCInBoth   is now obsolete!");
        exit(1);
        // const Float_t kZbegADC1 = -kZbegFrontBar-2.;  // (ecalvovi@cern.ch) 
        // const Float_t kZbegADC2 = -1885.0;            // (ecalvovi@cern.ch) 
        // ad -> AddNode(vADCarray , 1, new TGeoTranslation(0., 0., kZbegADC1 - kADCCelldz/2.)); // Tunnel
        // ad -> AddNode(vADCarrayH, 2, new TGeoTranslation(0., 0., kZbegADC2 - kADCCelldz/2.)); // Cavern 
        break;
      }
  }


  // ==========================================================================
  // 
  // Add structure volumes to top volume
  //
  // ==========================================================================

  TGeoVolumeAssembly * top = new TGeoVolumeAssembly("voADStruct");
  top->AddNode(voSaa3EndPlate, 1, new TGeoTranslation( 0., 0., kZendAbs + 1.95/2.));
  z = kZwall;
  top->AddNode(voWallBigPlate, 1, new TGeoTranslation(0., 0., z - 0.5 * dAlWallThick ));
  top->AddNode(voWallSqrPlate, 1, new TGeoTranslation(0., 0., z + 0.5 * dAlWallThick ));
  z = kZendAbs + 1.95 + dzRodL/2.; 
  top->AddNode(voSaa3Rod,  1, new TGeoTranslation(  12.5, -12.75, z));
  top->AddNode(voSaa3Rod,  2, new TGeoTranslation(  12.5,  12.75, z));
  top->AddNode(voSaa3Rod,  3, new TGeoTranslation( -12.5, -12.75, z));
  top->AddNode(voSaa3Rod,  4, new TGeoTranslation( -12.5,  12.75, z));
  top->AddNode(voValve,    1, new TGeoTranslation( 0., 0., zPosValve));
  //
  top->AddNode(voMoVMAOI,  1, new TGeoTranslation( 0., 0., kZbegVMAOI));
  //
  top->AddNode(voFrontBar, 1, new TGeoTranslation( 0., 0., kZbegFrontBar + dzF/2.));
  z = kZbegCoil;
  top->AddNode(voCoil, 1, new TGeoCombiTrans(  3.6 + dz/2., 0., z, Ry90m));
  top->AddNode(voCoil, 2, new TGeoCombiTrans(  3.6 + dz/2., 0., z, new TGeoRotation((*Ry90m)*(*Rx180))));
  top->AddNode(voCoil, 3, new TGeoCombiTrans( -3.6 - dz/2., 0., z, Ry90m));
  top->AddNode(voCoil, 4, new TGeoCombiTrans( -3.6 - dz/2., 0., z, new TGeoRotation((*Ry90m)*(*Rx180))));
  z = kZbegFrontBar + dzF + kdzLatBar/2.;
  top->AddNode(voLatBar, 1, new TGeoTranslation(  31.9, 0., z));
  top->AddNode(voLatBar, 2, new TGeoTranslation( -31.9, 0., z));
  //
  // Add structures (top) to AD node
  //
  if (fADCstruct) {
    ad->AddNode(top,1, Ry180);
  }

  //
  // Add Everything to ALICE
  //
  TGeoVolume *alice = gGeoManager->GetVolume("ALIC");
  alice->AddNode(ad, 1);
  

  // gGeoManager->DefaultColors();
  // gGeoManager->CloseGeometry();
  // gGeoManager->SetVisLevel(10);
  // gGeoManager->SetVisOption(0);
  // alice->Draw("ogl");

  return; 
  printf("<=== AliADv1::CreateAD(): ver=[Feb 3st, 2015]; contact=[ecalvovi@cern.ch]\n");
}

//_____________________________________________________________________________
void AliADv1::AddAlignableVolumes() const
{
   //
   // Create entries for alignable volumes associating the symbolic volume
   // name with the corresponding volume path. Needs to be syncronized with
   // eventual changes in the geometry.
   //
   // ADA and ADC 
   TString volpath1, volpath2, volpath3, volpath4;
   TString symname1, symname2, symname3, symname4;

   symname1 = "AD/ADC1";
   symname2 = "AD/ADC2"; 
   symname3 = "AD/ADA1";
   symname4 = "AD/ADA2"; 
   switch (fADCPosition) {
     case kADCInTunnel:
       volpath1 = "/ALIC_1/AD_1/ADC_1";
       volpath2 = "/ALIC_1/AD_1/ADC_2";
       break;
     case kADCInCavern:
       volpath1 = "/ALIC_1/AD_1/ADCh_1";
       volpath2 = "/ALIC_1/AD_1/ADCh_2";
       break;
     case kADCInBoth:
       volpath1 = "/ALIC_1/AD_1/ADC_1";
       volpath2 = "/ALIC_1/AD_1/ADCh_2";
       break;
   }
   volpath3 = "/ALIC_1/AD_1/ADA_1";
   volpath4 = "/ALIC_1/AD_1/ADA_2";

   if ( !gGeoManager->SetAlignableEntry(symname1.Data(), volpath1.Data()) )
      AliFatal(Form( "Alignable entry %s not created. Volume path %s not valid", symname1.Data(), volpath1.Data()) );
   if ( GetADCTwoInstalled() && !gGeoManager->SetAlignableEntry(symname2.Data(), volpath2.Data()) )
      AliFatal(Form( "Alignable entry %s not created. Volume path %s not valid", symname2.Data(), volpath2.Data()) );
   if ( !gGeoManager->SetAlignableEntry(symname3.Data(), volpath3.Data()) )
      AliFatal(Form( "Alignable entry %s not created. Volume path %s not valid", symname3.Data(), volpath3.Data()) );
   if ( GetADATwoInstalled() && !gGeoManager->SetAlignableEntry(symname4.Data(), volpath4.Data()) )
      AliFatal(Form("Alignable entry %s not created. Volume path %s not valid", symname4.Data(), volpath4.Data()) );
   
}


//_____________________________________________________________________________
void AliADv1::StepManager()
{

   //
   // Routine called at every step in the AD
   //

   // ADA and ADC static Variables         //
   // static  Int_t   numStep_ad = 0;      //
   // static  Int_t   vol_ad[2];           //
  
   /////////////////////////////////////////////////////////////////////////
   //                            ADA and ADC                              //
   /////////////////////////////////////////////////////////////////////////
      
      
   // Get sensitive volumes id (scintillator pads)
   static Int_t idADA  = gMC->VolId( "ADApad" );
   static Int_t idADC  = gMC->VolId( "ADCpad" );
   // static Int_t idADCh = gMC->VolId( "ADCpadH" );
   
   static Bool_t fOnlyOnce = kTRUE;
   if (fOnlyOnce) {
     //printf("  gMC->VolId(\"ADApad\" ) = %3d\n", idADA);
     //printf("  gMC->VolId(\"ADCpad\" ) = %3d\n", idADC);
     // printf("  gMC->VolId(\"ADCpadH\") = %3d\n", idADCh);
     fOnlyOnce = kFALSE;
   }

   // We keep only charged tracks : 
   // if ( !gMC->TrackCharge() || !gMC->IsTrackAlive() ) return;   
   // We keep charged and non-charged tracks : 
   if ( !gMC->IsTrackAlive() ) return;   
   
   Int_t copy;
   Int_t current_volid = gMC->CurrentVolID( copy );

   // check is the track is in a sensitive volume
   // if( current_volid != idADA && current_volid != idADC && current_volid != idADCh ) {
   if( current_volid != idADA && current_volid != idADC ) {
      return; // not in the sensitive volume 
   }
   
   // First read the position, otherwise weird results! //ecv
   Double_t s[3];
   Float_t  x[3];
   gMC->TrackPosition( s[0], s[1], s[2] );
   for ( Int_t j=0; j<3; j++ ) x[j] = s[j];
   
   // Get sector copy (1,2,3,4) ( 1 level up from pad )
   Int_t sect;
   gMC->CurrentVolOffID( 1, sect );
   
   // Get Detector copy (1,2) ( 2 levels up from pad )
   Int_t detc;
   gMC->CurrentVolOffID( 2, detc );

   // Set detector type: ADA or ADC
   Int_t ADlayer = (current_volid == idADC ) ? 0 : 2;

   //printf("CurVolID: %d | sect: %2d | detc: %2d\n", current_volid, sect, detc); 

   sect--;          //     sector within layer [0-3]
   detc--;          //     detector copy       [0-1]
   ADlayer += detc; //     global layer number [0-3]

   Int_t ADsector = ADlayer*4 + sect; // Global AD sector number [0-15]
   // Layer    Sector Number 
   // ADC 0  =   0- 3
   // ADC 1  =   4- 7
   // ADA 2  =   8-11
   // ADA 3  =  12-15
   //printf("\n ADsector: %2d | ADlayer: %2d | sect: %2d | x: %8.2f | y: %8.2f | z: %8.2f\n", ADsector, ADlayer, sect, x[0], x[1], x[2]); // Debug ECV
   
   Double_t lightYield_ad;
   Double_t photoCathodeEfficiency;
  
   if( ADlayer <2 )  {
      lightYield_ad          = fADCLightYield;
      photoCathodeEfficiency = fADCPhotoCathodeEfficiency;
   } else  {
      lightYield_ad          = fADALightYield;
      photoCathodeEfficiency = fADAPhotoCathodeEfficiency;
   }
      
   Float_t destep_ad = gMC->Edep();
   Float_t step_ad   = gMC->TrackStep();
   Int_t  nPhotonsInStep_ad = Int_t( destep_ad / (lightYield_ad * 1e-9) ); 
   nPhotonsInStep_ad = gRandom->Poisson( nPhotonsInStep_ad );
   
   static  Float_t eloss_ad    = 0.;
   static  Float_t tlength_ad  = 0.;   
   static  Int_t   nPhotons_ad = 0;      
   static  Float_t hits_ad[11];            
   static  Int_t   vol_ad[5];

   eloss_ad   += destep_ad;
   tlength_ad += step_ad;  
 
   if ( gMC->IsTrackEntering() ) { 
      nPhotons_ad = nPhotonsInStep_ad;
      Double_t p[4];
      gMC->TrackMomentum( p[0], p[1], p[2], p[3] );
      Float_t pt  = TMath::Sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] ); 
      TParticle *par = gAlice->GetMCApp()->Particle(gAlice->GetMCApp()->GetCurrentTrackNumber());
      Int_t imo = par->GetFirstMother();
      Int_t pdgMo = 0;
      if ( imo > 0 ) {
         TParticle * pmot = gAlice->GetMCApp()->Particle(imo);
         pdgMo = pmot->GetPdgCode();
      }

      // Set integer values
      vol_ad[0]  = par->GetStatusCode();    // secondary flag //ecv
      vol_ad[1]  = par->GetPdgCode();       // PDG code
      vol_ad[2]  = pdgMo;                   // PDG of the mother
      // Set float values
      hits_ad[0]  = x[0];     // X
      hits_ad[1]  = x[1];     // Y 
      hits_ad[2]  = x[2];     // Z       
      hits_ad[3]  = p[3];     // kinetic energy of the entering particle
      hits_ad[4]  = pt;       // Pt
      hits_ad[5]  = p[0];     // Px
      hits_ad[6]  = p[1];     // Py
      hits_ad[7]  = p[2];     // Pz
      hits_ad[8]  = 1.0e09*gMC->TrackTime(); // in ns!
  
      tlength_ad = 0.0;
      eloss_ad   = 0.0; 
      
      return; // without return, we count 2 times nPhotonsInStep_ad !!!???
   }
   
   nPhotons_ad += nPhotonsInStep_ad;

   if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared() ) {

      // Set integer values
      vol_ad[3]  = nPhotons_ad;
      vol_ad[4]  = ADsector;  // sector number (scintillator ID)
      // Set float values
      hits_ad[9]  = tlength_ad;    // track lenght inside ADC or ADA
      hits_ad[10] = eloss_ad;      // energy loss
      Int_t track; 
      if(fKeepHistory) track = gAlice->GetMCApp()->GetCurrentTrackNumber();
      else track = gAlice->GetMCApp()->GetPrimary( gAlice->GetMCApp()->GetCurrentTrackNumber() );
      AddHit( track, vol_ad, hits_ad ); // <-- this is in AliAD.cxx
      tlength_ad        = 0.0;
      eloss_ad          = 0.0; 
      nPhotons_ad       = 0;
   }
       
   //   Do we need track reference ????
   // if( gMC->IsTrackEntering() || gMC->IsTrackExiting() ) {
   //    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), 49);
   // }
}
//_________________________________________________________
void AliADv1::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
	TClonesArray &lhits = *fHits;
	new(lhits[fNhits++]) AliADhit(fIshunt,track,vol,hits);
}
//_________________________________________________________
// void AliADv1::AddDigits(Int_t* track, Int_t module, Float_t time)
// {
// 	TClonesArray &ldigits = *fDigits;
// 	new(ldigits[fNdigits++]) AliADdigit(track,module,time);
// }
//_________________________________________________________
void AliADv1::MakeBranch(Option_t *option)
{

	// Create branches in the current tree
	TString branchname(Form("%s",GetName()));
	AliDebug(2,Form("fBufferSize = %d",fBufferSize));
	const char *cH = strstr(option,"H");
	if (fHits && fLoader->TreeH() && cH)
	{
		fLoader->TreeH()->Branch(branchname.Data(),&fHits,fBufferSize);
		AliDebug(2,Form("Making Branch %s for hits",branchname.Data()));
	}
	const char *cD = strstr(option,"D");
  	if (fDigits   && fLoader->TreeD() && cD) 
	{
    		fLoader->TreeD()->Branch(branchname.Data(),&fDigits, fBufferSize);
    		AliDebug(2,Form("Making Branch %s for digits",branchname.Data()));
  	}  
}
