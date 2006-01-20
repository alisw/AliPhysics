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
  */
  //End_Html
  
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
  //  Define Mother 
  //  
  TGeoPgon* shMother = new TGeoPgon(22.5, 360., 8, 2);
  shMother->DefineSection(0, -600., 580., 790.);
  shMother->DefineSection(1,  600., 580., 790.);  
  // 
  TGeoVolume* voMother = new TGeoVolume("L3MO", shMother, medAir);
  //  
  // Define coils 
  //
  TGeoPgon* shCoils = new TGeoPgon(22.5, 360., 8, 2);
  shCoils->DefineSection(0, -600., 585., 690.);
  shCoils->DefineSection(1,  600., 585., 690.);  
  // 
  TGeoVolume* voCoils = new TGeoVolume("L3C0", shCoils, medAlu);
  voMother->AddNode(voCoils, 1, new TGeoTranslation(0., 0., 0.));
  //
  TGeoPgon* shCoilsI = new TGeoPgon(22.5, 360., 8, 2);
  shCoilsI->DefineSection(0, -600., 580., 585.);
  shCoilsI->DefineSection(1,  600., 580., 585.);  
  //
  TGeoVolume* voCoilsI = new TGeoVolume("L3C1", shCoilsI, medAluI);
  voMother->AddNode(voCoilsI, 1, new TGeoTranslation(0., 0., 0.));
  //
  // Define yoke
  //
  TGeoPgon* shYoke = new TGeoPgon(22.5, 360., 8, 2);
  shYoke->DefineSection(0, -600., 690., 790.);
  shYoke->DefineSection(1,  600., 690., 790.);  
  // 
  TGeoVolume* voYoke = new TGeoVolume("L3YO", shYoke, medFe);
  voMother->AddNode(voYoke, 1, new TGeoTranslation(0., 0., 0.));

  //
  // Define the return yoke of L3 ("Doors") 
  //
  // Original outer part
  TGeoPgon* shDoorO = new TGeoPgon(22.5, 360., 8, 2);
  shDoorO->DefineSection(0,  600., 240., 790.);
  shDoorO->DefineSection(1,  700., 240., 790.);  
  shDoorO->SetName("A");
  //
  // Additional inner part
  TGeoPgon* shDoorI = new TGeoPgon(22.5, 360., 8, 3);
  shDoorI->DefineSection(0,  600., 163.5, 270.);
  shDoorI->DefineSection(1,  670., 163.5, 270.);  
  shDoorI->DefineSection(2,  700., 213.5, 270.);
  shDoorI->SetName("B"); 
  //
  // For transport: low thresholds close to chambers requires special volume
  //
  TGeoPgon* shDoorIe = new TGeoPgon(22.5, 360., 8, 3);
  shDoorIe->DefineSection(0,  600., 163.5, 168.5);
  shDoorIe->DefineSection(1,  670., 163.5, 168.5);  
  shDoorIe->DefineSection(2,  700., 213.5, 218.5);
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
  //
  // The assembly of everything
  //
  TGeoVolumeAssembly* l3 = new TGeoVolumeAssembly("L3TV");
  TGeoRotation* rotxz = new TGeoRotation("rotxz",  90., 0., 90., 90., 180., 0.);
  l3->AddNode(voMother, 1, new TGeoTranslation(0., 0., 0.));
  l3->AddNode(voDoor, 1, new TGeoTranslation(0., 0., 0.));  
  l3->AddNode(voDoor, 2, new TGeoCombiTrans(0., 0., 0., rotxz));
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

