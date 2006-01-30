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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  L3 Magnet                                                                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliMAGClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:andreas.morsch@cern.ch">Andreas Morsch</a>.
</font>
<pre>

*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include <TVirtualMC.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TGeoPgon.h>
#include <TGeoCompositeShape.h>
#include <TGeoManager.h>

#include "AliMAG.h"
#include "AliMagF.h"
#include "AliRun.h"
 
ClassImp(AliMAG)
 
//_____________________________________________________________________________
AliMAG::AliMAG()
{
  //
  // Default constructor for L3 magnet
  //
}
 
//_____________________________________________________________________________
AliMAG::AliMAG(const char *name, const char *title)
  : AliModule(name,title)
{
  //
  // Standard constructor for L3 magnet
  //
  //Begin_Html
  /*
    <img src="picts/aliMAG.gif">
  */
  //End_Html
  
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}

//_____________________________________________________________________________
void AliMAG::CreateGeometry()
{
  //
  // Create geometry for L3 magnet
  //
  //Begin_Html
  /*
    <img src="picts/mag.gif">
  */
  //End_Html
    
  //Begin_Html
  /*
    <img src="picts/tree_mag.gif">
    <br> Dimensions taken from drawing: ALIL3___00010
  //End_Html
  */
// Octagon
    const Int_t   kNSides              =    8;
    const Float_t kStartAngle          =   22.5; // deg
    const Float_t kFullAngle           =  360.0; // deg
//  Mother volume 
    const Float_t kRBMotherInner       = 560.00; // cm
    const Float_t kRBMotherOuter       = 790.50; // cm
    const Float_t kLBMother            = 706.00; // cm
// Yoke     
    const Float_t kRYokeInner          = 703.50; // cm
    const Float_t kRYokeOuter          = 790.50; // cm
    const Float_t kLYoke               = 620.00; // cm
// Coil
    const Float_t kRCoilInner          = 593.00; // cm
    const Float_t kRCoilOuter          = 682.00; // cm
    const Float_t kLCoil               = 587.30; // cm
// Thermal Shield    
    const Float_t kRThermalShieldInner = 566.00; // cm
    const Float_t kRThermalShieldOuter = 571.00; // cm
// Crown    
    const Float_t kRCrownInner         = 560.00; // cm    
    const Float_t kRCrownOuter         = 785.50; // cm
    const Float_t kLCrown1             = 605.00; // cm
    const Float_t kLCrown2             = 620.00; // cm
    const Float_t kLCrown3             = 706.00; // cm
// Door
    const Float_t kRDoorInner          = 246.50; // cm
    const Float_t kRDoorOuter          = 560.00; // cm
    const Float_t kLDoor1              = 615.50; // cm
    const Float_t kLDoor2              = 714.60; // cm
    
    
  //
  // Top volume 
  TGeoVolume* top = gGeoManager->GetVolume("ALIC");
  // Media 
  TGeoMedium* medAir  = gGeoManager->GetMedium("MAG_AIR_C1");
  TGeoMedium* medAlu  = gGeoManager->GetMedium("MAG_ALU_C1");  
  TGeoMedium* medAluI = gGeoManager->GetMedium("MAG_ALU_C0");
  TGeoMedium* medFe   = gGeoManager->GetMedium("MAG_FE_C1");
  TGeoMedium* medFeI  = gGeoManager->GetMedium("MAG_FE_C0");
  //
  // Offset between LHC and LEP axis
  Float_t os = -30.;

  //
  //  Define Barrel Mother 
  //  
  TGeoPgon* shBMother = new TGeoPgon(kStartAngle, kFullAngle, kNSides, 2);
  shBMother->DefineSection(0, -kLBMother, kRBMotherInner, kRBMotherOuter);
  shBMother->DefineSection(1,  kLBMother, kRBMotherInner, kRBMotherOuter);  
  // 
  TGeoVolume* voBMother = new TGeoVolume("L3BM", shBMother, medAir);
  //
  // Define Thermal Shield
  //
  TGeoPgon* shThermSh = new TGeoPgon(kStartAngle, kFullAngle, kNSides, 2);
  shThermSh->DefineSection(0, -kLCoil, kRThermalShieldInner, kRThermalShieldOuter);
  shThermSh->DefineSection(1,  kLCoil, kRThermalShieldInner, kRThermalShieldOuter);  
  // 
  TGeoVolume* voThermSh = new TGeoVolume("L3TS", shThermSh, medAluI);
  voBMother->AddNode(voThermSh, 1, new TGeoTranslation(0., 0., 0.));
  //  
  // Define Coils
  //
  TGeoPgon* shCoils = new TGeoPgon(kStartAngle, kFullAngle, kNSides, 2);
  shCoils->DefineSection(0, -kLCoil, kRCoilInner, kRCoilOuter);
  shCoils->DefineSection(1,  kLCoil, kRCoilInner, kRCoilOuter);  
  // 
  TGeoVolume* voCoils = new TGeoVolume("L3C0", shCoils, medAlu);
  voBMother->AddNode(voCoils, 1, new TGeoTranslation(0., 0., 0.));
  //
  // Define Yoke
  //
  TGeoPgon* shYoke = new TGeoPgon(kStartAngle, kFullAngle, kNSides, 2);
  shYoke->DefineSection(0, -kLYoke, kRYokeInner, kRYokeOuter);
  shYoke->DefineSection(1, +kLYoke, kRYokeInner, kRYokeOuter);  
  // 
  TGeoVolume* voYoke = new TGeoVolume("L3YO", shYoke, medFe);
  voBMother->AddNode(voYoke, 1, new TGeoTranslation(0., 0., 0.));

  //
  // Define Crown
  //
  TGeoPgon* shCrown = new TGeoPgon(kStartAngle, kFullAngle, kNSides, 4);
  shCrown->DefineSection(0,  kLCrown1, kRCrownInner, kRYokeInner);
  shCrown->DefineSection(1,  kLCrown2, kRCrownInner, kRYokeInner);  
  shCrown->DefineSection(2,  kLCrown2, kRCrownInner, kRCrownOuter);  
  shCrown->DefineSection(3,  kLCrown3, kRCrownInner, kRCrownOuter);  
  // 
  TGeoVolume* voCrown = new TGeoVolume("L3CR", shCrown, medFe);
  //
  // Define Door 
  //
  // Original outer part
  TGeoPgon* shDoorO = new TGeoPgon(kStartAngle, kFullAngle, kNSides, 2);
  shDoorO->DefineSection(0,  kLDoor1, kRDoorInner, kRDoorOuter);
  shDoorO->DefineSection(1,  kLDoor2, kRDoorInner, kRDoorOuter);  
  shDoorO->SetName("A");
  //
  // Additional inner part
  TGeoPgon* shDoorI = new TGeoPgon(kStartAngle, kFullAngle, kNSides, 3);
  shDoorI->DefineSection(0,  kLDoor1, 163.5, 280.);
  shDoorI->DefineSection(1,  686.,    163.5, 280.);  
  shDoorI->DefineSection(2,  kLDoor2, 213.5, 280.);
  shDoorI->SetName("B"); 
  //
  // For transport: low thresholds close to chambers requires special medium
  //
  TGeoPgon* shDoorIe = new TGeoPgon(kStartAngle, kFullAngle, kNSides, 3);
  shDoorIe->DefineSection(0,  kLDoor1, 163.5, 168.5);
  shDoorIe->DefineSection(1,  686.,    163.5, 168.5);  
  shDoorIe->DefineSection(2,  kLDoor2, 213.5, 218.5);
  TGeoVolume* voDoorIe = new TGeoVolume("L3DE", shDoorIe, medFeI);
  //
  // Use composite shape here to account for the excentric door opening.
  // This avoids the overlap with the beam shield and the muon tracking station 1
  //
  TGeoTranslation* offset = new TGeoTranslation("t1", 0., -os, 0.);
  offset->RegisterYourself();
    
  TGeoCompositeShape* shDoor = new TGeoCompositeShape("L3Door", "A+B:t1");
  //
  TGeoVolume* voDoor = new TGeoVolume("L3DO", shDoor, medFe);
  voDoor->AddNode(voDoorIe, 1, new TGeoTranslation(0., -os, 0.));
  // Position crown and door
  TGeoRotation* rotxz = new TGeoRotation("rotxz",  90., 0., 90., 90., 180., 0.);

  TGeoVolumeAssembly *l3 = new TGeoVolumeAssembly("L3MO");
  voBMother->AddNode(voCrown, 1, new TGeoTranslation(0., 0., 0.));  
  voBMother->AddNode(voCrown, 2, new TGeoCombiTrans(0., 0., 0., rotxz));
  l3->AddNode(voBMother, 1, new TGeoTranslation(0.,0.,0.));
  l3->AddNode(voDoor,  1, new TGeoTranslation(0., 0., 0.));  
  l3->AddNode(voDoor,  2, new TGeoCombiTrans(0., 0., 0., rotxz));
  top->AddNode(l3, 1, new TGeoTranslation(0., os, 0.));
}

//_____________________________________________________________________________
void AliMAG::CreateMaterials()
{
  //
  // Create materials for L3 magnet
  //
  
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  Float_t epsil, stmin, deemax, tmaxfd, stemax;


  // --- Define the various materials for GEANT --- 
  
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;


  //     Aluminum 
  AliMaterial(9, "Al0$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(29, "Al1$", 26.98, 13., 2.7, 8.9, 37.2);
  
  //     Iron 
  AliMaterial(10, "Fe0$", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(30, "Fe1$", 55.85, 26., 7.87, 1.76, 17.1);
  
  //     Air 
  AliMixture(15, "AIR0$      ", aAir, zAir, dAir, 4, wAir);
  AliMixture(35, "AIR1$      ", aAir, zAir, dAir, 4, wAir);
  
  // **************** 
  //     Defines tracking media parameters. 
  //     Les valeurs sont commentees pour laisser le defaut 
  //     a GEANT (version 3-21, page CONS200), f.m. 
  epsil  = .001;  // Tracking precision, 
  stemax = -1.;   // Maximum displacement for multiple scat 
  tmaxfd = -20.;  // Maximum angle due to field deflection 
  deemax = -.3;   // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  // *************** 
  
  //    IRON 
  
  AliMedium(10, "FE_C0             ", 10, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(30, "FE_C1             ", 30, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //     ALUMINUM 

  AliMedium(9, "ALU_C0            ",  9, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(29, "ALU_C1            ", 29, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //     AIR 
  
  AliMedium(15, "AIR_C0            ", 15, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(35, "AIR_C1            ", 35, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
}

//_____________________________________________________________________________
void AliMAG::DrawModule() const
{
  //
  // Draw a shaded view of the L3 magnet
  //
}

//_____________________________________________________________________________
void AliMAG::Init()
{
  //
  // Initialise L3 magnet after it has been built
  Int_t i;
  //
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" MAG_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    // Here the MAG initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}

