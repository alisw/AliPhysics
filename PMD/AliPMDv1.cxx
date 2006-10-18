/***************************************************************************
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

//
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Photon Multiplicity Detector Version 1                                   //
//  Bedanga Mohanty : February 14th 2006
//                                                                           //
//Begin_Html
/*
<img src="picts/AliPMDv1Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
////

#include "Riostream.h"
#include <TVirtualMC.h>

#include "AliConst.h" 
#include "AliLog.h"
#include "AliMC.h"
#include "AliMagF.h" 
#include "AliPMDv1.h"
#include "AliRun.h"

const Int_t   AliPMDv1::fgkNcolUM1    = 48;  // Number of cols in UM, type 1
const Int_t   AliPMDv1::fgkNcolUM2    = 96;  // Number of cols in UM, type 2
const Int_t   AliPMDv1::fgkNrowUM1    = 96;  // Number of rows in UM, type 1
const Int_t   AliPMDv1::fgkNrowUM2    = 48;  // Number of rows in UM, type 2
const Float_t AliPMDv1::fgkCellRadius = 0.25;     // Radius of a hexagonal cell
const Float_t AliPMDv1::fgkCellWall   = 0.02;     // Thickness of cell Wall
const Float_t AliPMDv1::fgkCellDepth  = 0.50;     // Gas thickness
const Float_t AliPMDv1::fgkThBase     = 0.2;      // Thickness of Base plate
const Float_t AliPMDv1::fgkThBKP      = 0.1;      // Thickness of Back plane
const Float_t AliPMDv1::fgkThAir      = 1.03;      // Thickness of Air
const Float_t AliPMDv1::fgkThPCB      = 0.16;     // Thickness of PCB
const Float_t AliPMDv1::fgkThLead     = 1.5;      // Thickness of Pb
const Float_t AliPMDv1::fgkThSteel    = 0.5;      // Thickness of Steel
const Float_t AliPMDv1::fgkGap        = 0.025;    // Air Gap
const Float_t AliPMDv1::fgkZdist      = 361.5;    // z-position of the detector
const Float_t AliPMDv1::fgkSqroot3    = 1.7320508;// Square Root of 3
const Float_t AliPMDv1::fgkSqroot3by2 = 0.8660254;// Square Root of 3 by 2
const Float_t AliPMDv1::fgkSSBoundary = 0.3;
const Float_t AliPMDv1::fgkThSS       = 1.03;
const Float_t AliPMDv1::fgkThG10      = 1.03;
ClassImp(AliPMDv1)
 
//_____________________________________________________________________________
AliPMDv1::AliPMDv1():
  fSMthick(0.),
  fDthick(0.),
  fSMLengthax(0.),
  fSMLengthay(0.),
  fSMLengthbx(0.),
  fSMLengthby(0.),
  fMedSens(0)
{
  //
  // Default constructor 
  //
  for (Int_t i = 0; i < 3; i++)
    {
      fDboxmm1[i]  = 0.;
      fDboxmm12[i] = 0.;
      fDboxmm2[i]  = 0.;
      fDboxmm22[i] = 0.;
    }
}
 
//_____________________________________________________________________________
AliPMDv1::AliPMDv1(const char *name, const char *title):
  AliPMD(name,title),
  fSMthick(0.),
  fDthick(0.),
  fSMLengthax(0.),
  fSMLengthay(0.),
  fSMLengthbx(0.),
  fSMLengthby(0.),
  fMedSens(0)
{
  //
  // Standard constructor
  //
  for (Int_t i = 0; i < 3; i++)
    {
      fDboxmm1[i]  = 0.;
      fDboxmm12[i] = 0.;
      fDboxmm2[i]  = 0.;
      fDboxmm22[i] = 0.;
    }
}

//_____________________________________________________________________________
void AliPMDv1::CreateGeometry()
{
  // Create geometry for Photon Multiplicity Detector

  GetParameters();
  CreateSupermodule();
  CreatePMD();
}

//_____________________________________________________________________________
void AliPMDv1::CreateSupermodule()
{
  // 
  // Creates the geometry of the cells of PMD, places them in  supermodule 
  // which is a rectangular object.
  // Basic unit is ECAR, a hexagonal cell made of Ar+CO2, which is 
  // placed inside another hexagonal cell made of Cu (ECCU) with larger 
  // radius, compared to ECAR. The difference in radius gives the dimension 
  // of half width of each cell wall.
  // These cells are placed in a rectangular strip which are of 2 types 
  // EST1 and EST2 
  // 2 types of unit modules are made EUM1 and EUM2 which contains these strips
  // placed repeatedly 
  // Each supermodule (ESMA, ESMB), made of G10 is filled with following 
  //components. They have 6 unit moudles inside them
  // ESMA, ESMB are placed in EPMD along with EMPB (Pb converter) 
  // and EMFE (iron support) 

  
  Int_t i,j;
  Int_t number;
  Int_t ihrotm,irotdm;
  Float_t xb, yb, zb;

  Int_t *idtmed = fIdtmed->GetArray()-599;
 
  AliMatrix(ihrotm, 90., 30.,   90.,  120., 0., 0.);
  AliMatrix(irotdm, 90., 180.,  90.,  270., 180., 0.);
 
  // STEP - I
  //******************************************************//
  // First create the sensitive medium of a hexagon cell (ECAR)
  // Inner hexagon filled with gas (Ar+CO2)
  
  Float_t hexd2[10] = {0.,360.,6,2,-0.25,0.,0.23,0.25,0.,0.23};
  hexd2[4] = -fgkCellDepth/2.;
  hexd2[7] =  fgkCellDepth/2.;
  hexd2[6] =  fgkCellRadius - fgkCellWall;
  hexd2[9] =  fgkCellRadius - fgkCellWall;
  
  gMC->Gsvolu("ECAR", "PGON", idtmed[604], hexd2,10);
  gMC->Gsatt("ECAR", "SEEN", 0);
  //******************************************************//

  // STEP - II
  //******************************************************//
  // Place the sensitive medium inside a hexagon copper cell (ECCU)
  // Outer hexagon made of Copper
  
  Float_t hexd1[10] = {0.,360.,6,2,-0.25,0.,0.25,0.25,0.,0.25};
  hexd1[4] = -fgkCellDepth/2.;
  hexd1[7] =  fgkCellDepth/2.;
  hexd1[6] =  fgkCellRadius;
  hexd1[9] =  fgkCellRadius;

  gMC->Gsvolu("ECCU", "PGON", idtmed[614], hexd1,10);
  gMC->Gsatt("ECCU", "SEEN", 0);
  gMC->Gsatt("ECCU", "COLO", 4);

  // Place  inner hex (sensitive volume) inside outer hex (copper)
  
  gMC->Gspos("ECAR", 1, "ECCU", 0., 0., 0., 0, "ONLY");
  //******************************************************//

  // STEP - III
  //******************************************************//
  // Now create Rectangular TWO strips (EST1, EST2) 
  // of 1 column and 48 or 96 cells length

  // volume for first strip EST1 made of AIR 

  Float_t dbox1[3];
  dbox1[0] = fgkCellRadius/fgkSqroot3by2;
  dbox1[1] = fgkNrowUM1*fgkCellRadius;
  dbox1[2] = fgkCellDepth/2.;
  
  gMC->Gsvolu("EST1","BOX", idtmed[698], dbox1, 3);
  gMC->Gsatt("EST1", "SEEN", 0);

  // volume for second strip EST2 


  Float_t dbox2[3];
  dbox2[1] = fgkNrowUM2*fgkCellRadius;
  dbox2[0] = dbox1[0];
  dbox2[2] = dbox1[2];

  gMC->Gsvolu("EST2","BOX", idtmed[698], dbox2, 3);
  gMC->Gsatt("EST2", "SEEN", 0);

  // Place hexagonal cells ECCU placed inside EST1 
  xb = 0.; 
  zb = 0.;
  yb = (dbox1[1]) - fgkCellRadius; 
  for (i = 1; i <= fgkNrowUM1; ++i) 
    {
      number = i;
      gMC->Gspos("ECCU", number, "EST1", xb,yb,zb, 0, "ONLY");
      yb -= (fgkCellRadius*2.);
    }

  // Place hexagonal cells ECCU placed inside EST2 
  xb = 0.; 
  zb = 0.;
  yb = (dbox2[1]) - fgkCellRadius; 
  for (i = 1; i <= fgkNrowUM2; ++i) 
    {
      number = i;
      gMC->Gspos("ECCU", number, "EST2", xb,yb,zb, 0, "ONLY");
      //PH      cout << "ECCU in EST2 ==> " << number << "\t"<<xb <<  "\t"<<yb <<endl;
      yb -= (fgkCellRadius*2.);
    }


  //******************************************************//
 

  // STEP - IV
  //******************************************************//
 // 2 types of rectangular shaped unit modules EUM1 and EUM2 (defined by BOX) 
  //---------------------------------EHC1 Start----------------------//
  // Create EHC1 : The honey combs for a unit module type 1
  // First step is to create a honey comb unit module.
  // This is named as EHC1, we will lay the EST1 strips of
  // honey comb cells inside it.
  
  //Dimensions of EHC1
  //X-dimension = Number of columns + cell radius
  //Y-dimension = Number of rows * cell radius/sqrt3by2 - (some factor)
  //Z-dimension = cell depth/2

  Float_t dbox3[3];
  dbox3[0] = (dbox1[0]*fgkNcolUM1)-(fgkCellRadius*fgkSqroot3*(fgkNcolUM1-1)/6.);   
  dbox3[1] = dbox1[1]+fgkCellRadius/2.;
  dbox3[2] = fgkCellDepth/2.;

  //Create a BOX, Material AIR
  gMC->Gsvolu("EHC1","BOX", idtmed[698], dbox3, 3);
  gMC->Gsatt("EHC1", "SEEN", 0);  
  // Place rectangular strips EST1 inside EHC1 unit module
  xb = dbox3[0]-dbox1[0];  
  
  for (j = 1; j <= fgkNcolUM1; ++j)  
    {
      if(j%2 == 0)
	{
	  yb = -fgkCellRadius/2.0;
	}
      else
	{
	  yb = fgkCellRadius/2.0;
	}
      number = j;
      gMC->Gspos("EST1",number, "EHC1", xb, yb , 0. , 0, "MANY");
      //The strips are being placed from top towards bottom of the module
      //This is because the first cell in a module in hardware is the top
      //left corner cell
      xb = (dbox3[0]-dbox1[0])-j*fgkCellRadius*fgkSqroot3;

    }
  //--------------------EHC1 done----------------------------------//


  //---------------------------------EHC2 Start----------------------//
  // Create EHC2 : The honey combs for a unit module type 2
  // First step is to create a honey comb unit module.
  // This is named as EHC2, we will lay the EST2 strips of
  // honey comb cells inside it.

  //Dimensions of EHC2
  //X-dimension = Number of columns + cell radius
  //Y-dimension = Number of rows * cell radius/sqrt3by2 - (some factor)
  //Z-dimension = cell depth/2

  dbox3[0] = (dbox1[0]*fgkNcolUM1)-(fgkCellRadius*fgkSqroot3*(fgkNcolUM1-1)/6.);   
  dbox3[1] = dbox1[1]+fgkCellRadius/2.;
  dbox3[2] = fgkCellDepth/2.;

  Float_t dbox4[3];

  dbox4[0] =(dbox2[0]*fgkNcolUM2)-(fgkCellRadius*fgkSqroot3*(fgkNcolUM2-1)/6.); 
  dbox4[1] = dbox2[1] + fgkCellRadius/2.;
  dbox4[2] = dbox3[2];
  
  //Create a BOX of AIR
  gMC->Gsvolu("EHC2","BOX", idtmed[698], dbox4, 3);
  gMC->Gsatt("EHC2", "SEEN", 0);

  // Place rectangular strips EST2 inside EHC2 unit module
  xb = dbox4[0]-dbox2[0]; 
  for (j = 1; j <= fgkNcolUM2; ++j) 
  {
    if(j%2 == 0)
  {
    yb = -fgkCellRadius/2.0;
  }
    else
  {
    yb = +fgkCellRadius/2.0;
  }
    number = j;
    gMC->Gspos("EST2",number, "EHC2", xb, yb , 0. ,0, "MANY");
    xb = (dbox4[0]-dbox2[0])-j*fgkCellRadius*fgkSqroot3;
  }
  

  //--------------------EHC2 done----------------------------------//


  // Now the job is to assmeble an Unit module
  // It will have the following components
  // (a) Base plate of G10 of 0.2 cm 
  // (b) Air gap  of 0.05 cm 
  // (c) Bottom PCB of 0.16 cm G10
  // (d) Honey comb 0f 0.5 cm
  // (e) Top PCB  of 0.16 cm G10
  // (f) Air gap of 0.16 cm
  // (g) Back Plane of 0.1 cm G10
  // (h) Then all around then we have an air gap of 0.5mm
  // (i) Then all around 0.5mm thick G10 insulation
  // (h) Then all around Stainless Steel boundary channel 0.3 cm thick
  //Let us first create them one by one
  //---------------------------------------------------//

  // ---------------- Lets do it first for UM Type A -----//

 //--------------------------------------------------//
  //Bottom and Top PCB : EPCA
  //===========================
  // Make a 1.6mm thick G10 Bottom and Top PCB for Unit module A
  // X-dimension same as EHC1 - dbox3[0]
  // Y-dimension same as EHC1 - dbox3[1]
  // Z-dimension 0.16/2 = 0.08 cm
  //-------------------------------------------------//
  Float_t dboxPcbA[3];
  dboxPcbA[0]      = dbox3[0]; 
  dboxPcbA[1]      = dbox3[1];       
  dboxPcbA[2]      = fgkThPCB/2.;
  
  //Top and Bottom PCB is a BOX of Material G10
  gMC->Gsvolu("EPCA","BOX", idtmed[607], dboxPcbA, 3);
  gMC->Gsatt("EPCA", "SEEN", 0);
  //--------------------------------------------------------//  
  //Back Plane : EBKA
  //==================
  // Make a 1.0mm thick Back Plane PCB for Unit module A
  // X-dimension same as EHC1 - dbox3[0]
  // Y-dimension same as EHC1 - dbox3[1]
  // Z-dimension 0.1/2 = 0.05 cm
  //------------------------------------------------------//
  Float_t dboxBPlaneA[3];
  dboxBPlaneA[0]   = dbox3[0]; 
  dboxBPlaneA[1]   = dbox3[1];       
  dboxBPlaneA[2]   = fgkThBKP/2.;
  
  //Back PLane PCB of MAterial G10
  gMC->Gsvolu("EBKA","BOX", idtmed[607], dboxBPlaneA, 3);
  gMC->Gsatt("EBKA", "SEEN", 0);
  //-------------------------------------------------------------//  

 //---------- That was all in the Z -direction of Unit Module A----//

  //  Now lets us construct the boundary arround the Unit Module --//
  // This boundary has 
  // (a) 0.5 mm X and Y and 10.3 mm Z dimension  AIR gap
  // (b) 0.5 mm X and Y and 10.3 mm Z dimension G10
  // (c) 3.0 mm X and Y and 12.3 mm Z dimension Stainless Steel



  //-------------------------------------------------//
  //AIR GAP between UM and Boundary : ECGA FOR PRESHOWER PLANE
  //==========================================================
  // Make a 10.3mm thick Air gap for Unit module A
  // X-dimension same as EHC1+0.05
  // Y-dimension same as EHC1+0.05
  // Z-dimension 1.03/2 = 0.515 cm
  Float_t dboxAir3A[3];
  dboxAir3A[0]         = dbox3[0]+(2.0*fgkGap); 
  dboxAir3A[1]         = dbox3[1]+(2.0*fgkGap); 
  dboxAir3A[2]         = fgkThAir/2.;

  //FOR PRESHOWER
  //Air gap is a BOX of Material Air
  gMC->Gsvolu("ECGA","BOX", idtmed[698], dboxAir3A, 3);
  gMC->Gsatt("ECGA", "SEEN", 0);

  //FOR VETO
  //Air gap is a BOX of Material Air
  gMC->Gsvolu("ECVA","BOX", idtmed[698], dboxAir3A, 3);
  gMC->Gsatt("ECVA", "SEEN", 0);
  //-------------------------------------------------//  

 //-------------------------------------------------//
  //G10 boundary between honeycomb and SS : EDGA
  //================================================
  // Make a 10.3mm thick G10 Boundary for Unit module A
  // X-dimension same as EHC1+Airgap+0.05
  // Y-dimension same as EHC1+Airgap+0.05
  // Z-dimension 1.03/2 = 0.515 cm
  Float_t dboxGGA[3];
  dboxGGA[0]         = dboxAir3A[0]+(2.0*fgkGap); 
  dboxGGA[1]         = dboxAir3A[1]+(2.0*fgkGap); 
  dboxGGA[2]         = fgkThG10/2.;

  //FOR PRESHOWER
  //G10 BOX 
  gMC->Gsvolu("EDGA","BOX", idtmed[607], dboxGGA, 3);
  gMC->Gsatt("EDGA", "SEEN", 0);

  //FOR VETO
  //G10 BOX 
  gMC->Gsvolu("EDVA","BOX", idtmed[607], dboxGGA, 3);
  gMC->Gsatt("EDVA", "SEEN", 0);

  //-------------------------------------------------//  
  //----------------------------------------------------------//
  //Stainless Steel Bounadry : ESSA
  //==================================
  // Make a 10.3mm thick Stainless Steel boundary for Unit module A
  // X-dimension same as EHC1 + Airgap + G10 + 0.3
  // Y-dimension same as EHC1 + Airgap + G10 + 0.3
  // Z-dimension 1.03/2 = 0.515 cm
  //------------------------------------------------------//
  // A Stainless Steel Boundary Channel to house the unit module

  Float_t dboxSS1[3];
  dboxSS1[0]           = dboxGGA[0]+fgkSSBoundary; 
  dboxSS1[1]           = dboxGGA[1]+fgkSSBoundary;       
  dboxSS1[2]           = fgkThSS/2.;
  
  //FOR PRESHOWER

  //Stainless Steel boundary - Material Stainless Steel
  gMC->Gsvolu("ESSA","BOX", idtmed[618], dboxSS1, 3);
  gMC->Gsatt("ESSA", "SEEN", 0);

  //FOR VETO
  //Stainless Steel boundary - Material Stainless Steel
  gMC->Gsvolu("ESVA","BOX", idtmed[618], dboxSS1, 3);
  gMC->Gsatt("ESVA", "SEEN", 0);

  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  // Here we need to place the volume in order ESSA -> EDGA -> ECGA 
  // this makes the SS boundary and the 0.5mm thick FR4 insulation in place, 
  // and the air volume ECGA acts as mother for the rest of components.
  // The above placeemnt is done at (0.,0.,0.) relative coordiante 
  // Now we place bottom PCB, honeycomb, top PCB in this volume. We donot place
  // unnecessary air volumes now. Just leave the gap as we are placing them
  // in  air only. This also reduces the number of volumes for geant to track.

// Tree structure for different volumes
//
//				EUM1
//				 |
//			--------------------
//			|        |         |
//		      EBPA	ESSA	  EBKA
//				 |
//				EDGA
//				 |
//				ECGA
//				 |
//			--------------------
//			|        |	   |
//		      EPCA(1)   EHC1	 EPCA(2)
//		     (bottom)	 |	(top PCB)
//				 |
//			    Sensitive volume
//				(gas)
//	


  //FOR VETO
//Creating the side channels 
// SS boundary channel, followed by G10 and Air Gap  
  gMC->Gspos("EDVA", 1, "ESVA", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("ECVA", 1, "EDVA", 0., 0., 0., 0, "ONLY");

//FOR PRESHOWER
  gMC->Gspos("EDGA", 1, "ESSA", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("ECGA", 1, "EDGA", 0., 0., 0., 0, "ONLY");

 // now other components, using Bedanga's code, but changing the values.
  //Positioning Bottom PCB, Honey Comb abd Top PCB in AIR

  //For veto plane
  //Positioning the Bottom 0.16 cm PCB
  Float_t zbpcb = -dboxAir3A[2] + (2.0*fgkGap) + fgkThPCB/2.;
  gMC->Gspos("EPCA", 1, "ECVA", 0., 0., zbpcb, 0, "ONLY");
  //Positioning the Honey Comb 0.5 cm
  Float_t zhc = zbpcb + fgkThPCB/2. + fgkCellDepth/2.;
  gMC->Gspos("EHC1", 1, "ECVA", 0., 0., zhc, 0, "ONLY");
  //Positioning the Top PCB 0.16 cm
  Float_t ztpcb = zhc + fgkCellDepth/2 + fgkThPCB/2.;
  gMC->Gspos("EPCA", 2, "ECVA", 0., 0., ztpcb, 0, "ONLY");


  //For Preshower plane the ordering is reversed
  //Positioning the Bottom 0.16 cm PCB
  zbpcb = -dboxAir3A[2] + fgkThPCB + fgkThPCB/2.;
  gMC->Gspos("EPCA", 1, "ECGA", 0., 0., zbpcb, 0, "ONLY");
  //Positioning the Honey Comb 0.5 cm
  zhc = zbpcb + fgkThPCB/2. + fgkCellDepth/2.;
  gMC->Gspos("EHC1", 1, "ECGA", 0., 0., zhc, 0, "ONLY");
  //Positioning the Top PCB 0.16 cm
  ztpcb = zhc + fgkCellDepth/2 + fgkThPCB/2.;
  gMC->Gspos("EPCA", 2, "ECGA", 0., 0., ztpcb, 0, "ONLY");




 //--------------Now let us construct final UM ---------------//
  // We will do it as follows :
  // (i)  First make a UM of air. which will have dimensions
  //      of the SS boundary Channel (in x,y) and of height 13.3mm
  //(ii)  Then we will place all the components

  //----------------------------------------------------------//
  // A  unit module type A of Air
  // Dimensions of Unit Module same as SS boundary channel
  Float_t dboxUM1[3];
  dboxUM1[0] = dboxSS1[0];
  dboxUM1[1] = dboxSS1[1];
  dboxUM1[2] = fgkThSS/2. +0.15; // 0.15 added to accomodate Base Plate at
  // the bottom and the backplane PCB at the top.

  //FOR PRESHOWER
  //Create a Unit module of above dimensions Material : AIR
  gMC->Gsvolu("EUM1","BOX", idtmed[698], dboxUM1, 3);
  gMC->Gsatt("EUM1", "SEEN", 0);
  //FOR VETO
  gMC->Gsvolu("EUV1","BOX", idtmed[698], dboxUM1, 3);
  gMC->Gsatt("EUV1", "SEEN", 0);

  //----------------------------------------------------------------//

  //BASE PLATE : EBPA
  //==================
  // Make a 2mm thick G10 Base plate for Unit module A
  // Base plate is as big as the final UM dimensions that is as 
  // SS boundary channel
  Float_t dboxBaseA[3];
  dboxBaseA[0]       = dboxSS1[0];
  dboxBaseA[1]       = dboxSS1[1];       
  dboxBaseA[2]       = fgkThBase/2.;
  
  //Base Blate is a G10 BOX
  gMC->Gsvolu("EBPA","BOX", idtmed[607], dboxBaseA, 3);
  gMC->Gsatt("EBPA", "SEEN", 0);
  //----------------------------------------------------//  

  //FOR VETO
  //- Placing of all components of UM in AIR BOX EUM1--//
  //(1)   FIRST PUT THE BASE PLATE
  Float_t zbaseplate = -dboxUM1[2] + fgkThBase/2.;
  gMC->Gspos("EBPA", 1, "EUV1", 0., 0., zbaseplate, 0, "ONLY");

  //(2)   NEXT PLACING the SS BOX 
  Float_t zss = zbaseplate + fgkThBase/2. + fgkThSS/2.;
  gMC->Gspos("ESVA", 1, "EUV1", 0., 0., zss, 0, "ONLY");
  
  // (3) Positioning the Backplane PCB 0.1 cm
  Float_t zbkp = zss + fgkThSS/2. + fgkThBKP/2.;
  gMC->Gspos("EBKA", 1, "EUV1", 0., 0., zbkp, 0, "ONLY");

  //FOR PRESHOWER
  // (3) Positioning the Backplane PCB 0.1 cm
  zbkp = -dboxUM1[2] + fgkThBKP/2.;
  gMC->Gspos("EBKA", 1, "EUM1", 0., 0., zbkp, 0, "ONLY");

  //(2)   NEXT PLACING the SS BOX 
  zss = zbkp + fgkThBKP/2. + fgkThSS/2.;
  gMC->Gspos("ESSA", 1, "EUM1", 0., 0., zss, 0, "ONLY");
  
  //(1)   FIRST PUT THE BASE PLATE
  zbaseplate = zss + fgkThSS/2 + fgkThBase/2.;
  gMC->Gspos("EBPA", 1, "EUM1", 0., 0., zbaseplate, 0, "ONLY");
  //-------------------- UM Type A completed ------------------------//



  //-------------------- Lets do the same thing for UM type B -------//
 //--------------------------------------------------//
  //Bottom and Top PCB : EPCB
  //===========================
  // Make a 1.6mm thick G10 Bottom and Top PCB for Unit module B
  // X-dimension same as EHC2 - dbox4[0]
  // Y-dimension same as EHC2 - dbox4[1]
  // Z-dimension 0.16/2 = 0.08 cm
  //-------------------------------------------------//
  Float_t dboxPcbB[3];
  dboxPcbB[0]      = dbox4[0]; 
  dboxPcbB[1]      = dbox4[1];       
  dboxPcbB[2]      = fgkThPCB/2.;
  
  //Top and Bottom PCB is a BOX of Material G10
  gMC->Gsvolu("EPCB","BOX", idtmed[607], dboxPcbB, 3);
  gMC->Gsatt("EPCB", "SEEN", 0);
  //--------------------------------------------------------//  
  //Back Plane : EBKB
  //==================
  // Make a 1.0mm thick Back Plane PCB for Unit module B
  // X-dimension same as EHC2 - dbox4[0]
  // Y-dimension same as EHC2 - dbox4[1]
  // Z-dimension 0.1/2 = 0.05 cm
  //------------------------------------------------------//
  Float_t dboxBPlaneB[3];
  dboxBPlaneB[0]   = dbox4[0]; 
  dboxBPlaneB[1]   = dbox4[1];       
  dboxBPlaneB[2]   = fgkThBKP/2.;
  
  //Back PLane PCB of MAterial G10
  gMC->Gsvolu("EBKB","BOX", idtmed[607], dboxBPlaneB, 3);
  gMC->Gsatt("EBKB", "SEEN", 0);
  //-------------------------------------------------------------//  

 //---------- That was all in the Z -direction of Unit Module B----//

  //  Now lets us construct the boundary arround the Unit Module --//
  // This boundary has 
  // (a) 0.5 mm X and Y and 10.3 mm Z dimension  AIR gap
  // (b) 0.5 mm X and Y and 10.3 mm Z dimension G10
  // (c) 3.0 mm X and Y and 12.3 mm Z dimension Stainless Steel

  //-------------------------------------------------//
  //AIR GAP between UM and Boundary : ECGB
  //================================================
  // Make a 10.3mm thick Air gap for Unit module B
  // X-dimension same as EHC2+0.05
  // Y-dimension same as EHC2+0.05
  // Z-dimension 1.03/2 = 0.515 cm
  Float_t dboxAir3B[3];
  dboxAir3B[0]         = dbox4[0]+(2.0*fgkGap); 
  dboxAir3B[1]         = dbox4[1]+(2.0*fgkGap);       
  dboxAir3B[2]         = fgkThAir/2.;

  //PRESHOWER
  //Air gap is a BOX of Material Air
  gMC->Gsvolu("ECGB","BOX", idtmed[698], dboxAir3B, 3);
  gMC->Gsatt("ECGB", "SEEN", 0);
  //VETO
  gMC->Gsvolu("ECVB","BOX", idtmed[698], dboxAir3B, 3);
  gMC->Gsatt("ECVB", "SEEN", 0);

  //-------------------------------------------------//  

 //-------------------------------------------------//
  //G10 boundary between honeycomb and SS : EDGB
  //================================================
  // Make a 10.3mm thick G10 Boundary for Unit module B
  // X-dimension same as EHC2+Airgap+0.05
  // Y-dimension same as EHC2+Airgap+0.05
  // Z-dimension 1.03/2 = 0.515 cm
  Float_t dboxGGB[3];
  dboxGGB[0]         = dboxAir3B[0]+(2.0*fgkGap); 
  dboxGGB[1]         = dboxAir3B[1]+(2.0*fgkGap);      
  dboxGGB[2]         = fgkThG10/2.;

  //PRESHOWER
  //G10 BOX 
  gMC->Gsvolu("EDGB","BOX", idtmed[607], dboxGGB, 3);
  gMC->Gsatt("EDGB", "SEEN", 0);
  //VETO
  gMC->Gsvolu("EDVB","BOX", idtmed[607], dboxGGB, 3);
  gMC->Gsatt("EDVB", "SEEN", 0);
  //-------------------------------------------------//  
  //----------------------------------------------------------//
  //Stainless Steel Bounadry : ESSB
  //==================================
  // Make a 10.3mm thick Stainless Steel boundary for Unit module B
  // X-dimension same as EHC2 + Airgap + G10 + 0.3
  // Y-dimension same as EHC2 + Airgap + G10 + 0.3
  // Z-dimension 1.03/2 = 0.515 cm
  //------------------------------------------------------//
  // A Stainless Steel Boundary Channel to house the unit module

  Float_t dboxSS2[3];
  dboxSS2[0]  = dboxGGB[0] + fgkSSBoundary; 
  dboxSS2[1]  = dboxGGB[1] + fgkSSBoundary;       
  dboxSS2[2]  = fgkThSS/2.;
  
  //PRESHOWER
  //Stainless Steel boundary - Material Stainless Steel
  gMC->Gsvolu("ESSB","BOX", idtmed[618], dboxSS2, 3);
  gMC->Gsatt("ESSB", "SEEN", 0);
  //VETO
  gMC->Gsvolu("ESVB","BOX", idtmed[618], dboxSS2, 3);
  gMC->Gsatt("ESVB", "SEEN", 0);
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  // Here we need to place the volume in order ESSB -> EDGB -> ECGB 
  // this makes the SS boiundary and the 0.5mm thick FR4 insulation in place, 
  // and the air volume ECGB acts as mother for the rest of components.
  // The above placeemnt is done at (0.,0.,0.) relative coordiante 
  // Now we place bottom PCB, honeycomb, top PCB in this volume. We donot place
  // unnecessary air volumes now. Just leave the gap as we are placing them
  // in  air only. This also reduces the number of volumes for geant to track.

// Tree structure for different volumes
//
//				EUM2
//				 |
//			--------------------
//			|        |         |
//		      EBPB	ESSB	  EBKB
//				 |
//				EDGB
//				 |
//				ECGB
//				 |
//			--------------------
//			|        |	   |
//		      EPCB(1)   EHC2	 EPCB(2)
//		     (bottom)	 |	(top PCB)
//				 |
//			    Sensitive volume
//				(gas)
//	

//PRESHOWER
//Creating the side channels
// SS boundary channel, followed by G10 and Air Gap  
  gMC->Gspos("EDGB", 1, "ESSB", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("ECGB", 1, "EDGB", 0., 0., 0., 0, "ONLY");
  //VETO
  gMC->Gspos("EDVB", 1, "ESVB", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("ECVB", 1, "EDVB", 0., 0., 0., 0, "ONLY");

 // now other components, using Bedang's code, but changing the values.
  //Positioning Bottom PCB, Honey Comb abd Top PCB in AIR

  //VETO
  //Positioning the Bottom 0.16 cm PCB
  Float_t zbpcb2 = -dboxAir3B[2] + (2.0*fgkGap) + fgkThPCB/2.;
  gMC->Gspos("EPCB", 1, "ECVB", 0., 0., zbpcb2, 0, "ONLY");
  //Positioning the Honey Comb 0.5 cm
  Float_t zhc2 = zbpcb2 + fgkThPCB/2. + fgkCellDepth/2.;
  gMC->Gspos("EHC2", 1, "ECVB", 0., 0., zhc2, 0, "ONLY");
  //Positioning the Top PCB 0.16 cm
  Float_t ztpcb2 = zhc2 + fgkCellDepth/2 + fgkThPCB/2.;
  gMC->Gspos("EPCB", 2, "ECVB", 0., 0., ztpcb2, 0, "ONLY");

  //PRESHOWER
  //For preshower plane the ordering is reversed
  //Positioning the Bottom 0.16 cm PCB
  zbpcb2 = -dboxAir3B[2] + fgkThPCB + fgkThPCB/2.;
  gMC->Gspos("EPCB", 1, "ECGB", 0., 0., zbpcb2, 0, "ONLY");
  //Positioning the Honey Comb 0.5 cm
  zhc2 = zbpcb2 + fgkThPCB/2. + fgkCellDepth/2.;
  gMC->Gspos("EHC2", 1, "ECGB", 0., 0., zhc2, 0, "ONLY");
  //Positioning the Top PCB 0.16 cm
  ztpcb2 = zhc2 + fgkCellDepth/2 + fgkThPCB/2.;
  gMC->Gspos("EPCB", 2, "ECGB", 0., 0., ztpcb2, 0, "ONLY");



 //--------------Now let us construct final UM ---------------//
  // We will do it as follows :
  // (i)  First make a UM of air. which will have dimensions
  //      of the SS boundary Channel (in x,y) and of height 13.3mm
  //(ii)  Then we will place all the components

  //----------------------------------------------------------//
  // A  unit module type B of Air
  // Dimensions of Unit Module same as SS boundary channel

  Float_t dboxUM2[3];
  dboxUM2[0] = dboxSS2[0];
  dboxUM2[1] = dboxSS2[1];
  dboxUM2[2] = fgkThSS/2. +0.15; // 0.15 added to accomodate Base Plate at
  // the bottom and the backplane PCB at the top.

  //PRESHOWER
  //Create a Unit module of above dimensions Material : AIR
  gMC->Gsvolu("EUM2","BOX", idtmed[698], dboxUM2, 3);
  gMC->Gsatt("EUM2", "SEEN", 0);

  //VETO
  gMC->Gsvolu("EUV2","BOX", idtmed[698], dboxUM2, 3);
  gMC->Gsatt("EUV2", "SEEN", 0);
  //----------------------------------------------------------------//

  //BASE PLATE : EBPB
  //==================
  // Make a 2mm thick G10 Base plate for Unit module B
  // Base plate is as big as the final UM dimensions that is as 
  // SS boundary channel
  Float_t dboxBaseB[3];
  dboxBaseB[0]       = dboxSS2[0];
  dboxBaseB[1]       = dboxSS2[1];       
  dboxBaseB[2]       = fgkThBase/2.;
  
  //Base Blate is a G10 BOX
  gMC->Gsvolu("EBPB","BOX", idtmed[607], dboxBaseB, 3);
  gMC->Gsatt("EBPB", "SEEN", 0);
  //----------------------------------------------------//  

  //VETO
  //- Placing of all components of UM in AIR BOX EUM2--//
  //(1)   FIRST PUT THE BASE PLATE
  Float_t zbaseplate2 = -dboxUM2[2] + fgkThBase/2.;
  gMC->Gspos("EBPB", 1, "EUV2", 0., 0., zbaseplate2, 0, "ONLY");

  //(2)   NEXT PLACING the SS BOX 
  Float_t zss2 = zbaseplate2 + fgkThBase/2. + fgkThSS/2.;
  gMC->Gspos("ESVB", 1, "EUV2", 0., 0., zss2, 0, "ONLY");
  
  // (3) Positioning the Backplane PCB 0.1 cm
  Float_t zbkp2 = zss2 + fgkThSS/2. + fgkThBKP/2.;
  gMC->Gspos("EBKB", 1, "EUV2", 0., 0., zbkp2, 0, "ONLY");



  //FOR PRESHOWER
  // (3) Positioning the Backplane PCB 0.1 cm
  zbkp2 = -dboxUM2[2] + fgkThBKP/2.;
  gMC->Gspos("EBKB", 1, "EUM2", 0., 0., zbkp2, 0, "ONLY");

  //(2)   NEXT PLACING the SS BOX 
  zss2 = zbkp2 + fgkThBKP/2. + fgkThSS/2.;
  gMC->Gspos("ESSB", 1, "EUM2", 0., 0., zss2, 0, "ONLY");
  
  //(1)   FIRST PUT THE BASE PLATE
  zbaseplate2 = zss2 + fgkThSS/2 + fgkThBase/2.;
  gMC->Gspos("EBPB", 1, "EUM2", 0., 0., zbaseplate2, 0, "ONLY");
  //-------------------- UM Type B completed ------------------------//


  //--- Now we need to make Lead plates of UM dimension -----//

  /**************************/
  //----------------------------------------------------------//
  // The lead convertor is of unit module size
  // Dimensions of Unit Module same as SS boundary channel

  Float_t dboxPba[3];
  dboxPba[0] = dboxUM1[0];
  dboxPba[1] = dboxUM1[1];
  dboxPba[2] = fgkThLead/2.;
  // Lead of UM dimension
  gMC->Gsvolu("EPB1","BOX", idtmed[600], dboxPba, 3);
  gMC->Gsatt ("EPB1", "SEEN", 0);

  Float_t dboxPbb[3];
  dboxPbb[0] = dboxUM2[0];
  dboxPbb[1] = dboxUM2[1];
  dboxPbb[2] = fgkThLead/2.;
  // Lead of UM dimension
  gMC->Gsvolu("EPB2","BOX", idtmed[600], dboxPbb, 3);
  gMC->Gsatt ("EPB2", "SEEN", 0);

  //----------------------------------------------------------------//

  // 2 types of Rectangular shaped supermodules (BOX) 
  //each with 6 unit modules 
  
  // volume for SUPERMODULE ESMA 
  //Space added to provide a gapping for HV between UM's
  //There is a gap of 0.15 cm between two Modules (UMs)
  // in x-direction and 0.1cm along y-direction

  Float_t dboxSM1[3];
  dboxSM1[0] = 3.0*dboxUM1[0] + (2.0*0.075);
  dboxSM1[1] = 2.0*dboxUM1[1] +  0.05;
  dboxSM1[2] = dboxUM1[2];

  //FOR PRESHOWER  
  gMC->Gsvolu("ESMA","BOX", idtmed[698], dboxSM1, 3);
  gMC->Gsatt("ESMA", "SEEN", 0);
  
  //FOR VETO
  gMC->Gsvolu("EMVA","BOX", idtmed[698], dboxSM1, 3);
  gMC->Gsatt("EMVA", "SEEN", 0);

  //Position the 6 unit modules in EMSA
  Float_t xa1,xa2,xa3,ya1,ya2; 
  xa1 =  dboxSM1[0] - dboxUM1[0];
  xa2 = xa1 - dboxUM1[0] - 0.15 - dboxUM1[0];
  xa3 = xa2 - dboxUM1[0] - 0.15 - dboxUM1[0];
  ya1 = dboxSM1[1]  - dboxUM1[1];
  ya2 = ya1 - dboxUM1[1] - 0.1 - dboxUM1[1];

  //PRESHOWER
  gMC->Gspos("EUM1", 1, "ESMA", xa1, ya1, 0., 0, "ONLY");
  gMC->Gspos("EUM1", 2, "ESMA", xa2, ya1, 0., 0, "ONLY");
  gMC->Gspos("EUM1", 3, "ESMA", xa3, ya1, 0., 0, "ONLY");
  gMC->Gspos("EUM1", 4, "ESMA", xa1, ya2, 0., 0, "ONLY");
  gMC->Gspos("EUM1", 5, "ESMA", xa2, ya2, 0., 0, "ONLY");
  gMC->Gspos("EUM1", 6, "ESMA", xa3, ya2, 0., 0, "ONLY");

  //VETO
  gMC->Gspos("EUV1", 1, "EMVA", xa1, ya1, 0., 0, "ONLY");
  gMC->Gspos("EUV1", 2, "EMVA", xa2, ya1, 0., 0, "ONLY");
  gMC->Gspos("EUV1", 3, "EMVA", xa3, ya1, 0., 0, "ONLY");
  gMC->Gspos("EUV1", 4, "EMVA", xa1, ya2, 0., 0, "ONLY");
  gMC->Gspos("EUV1", 5, "EMVA", xa2, ya2, 0., 0, "ONLY");
  gMC->Gspos("EUV1", 6, "EMVA", xa3, ya2, 0., 0, "ONLY");


  // volume for SUPERMODULE ESMB 
  //Space is added to provide a gapping for HV between UM's
  Float_t dboxSM2[3];
  dboxSM2[0] = 2.0*dboxUM2[0] + 0.075; 
  dboxSM2[1] = 3.0*dboxUM2[1] + (2.0*0.05);
  dboxSM2[2] = dboxUM2[2];
  
  //PRESHOWER
  gMC->Gsvolu("ESMB","BOX", idtmed[698], dboxSM2, 3);
  gMC->Gsatt("ESMB", "SEEN", 0);
  //VETO 
  gMC->Gsvolu("EMVB","BOX", idtmed[698], dboxSM2, 3);
  gMC->Gsatt("EMVB", "SEEN", 0);

  //Position the 6 unit modules in EMSB
  Float_t xb1,xb2,yb1,yb2,yb3; 
  xb1 = dboxSM2[0] - dboxUM2[0];
  xb2 = xb1 - dboxUM2[0] - 0.15 - dboxUM2[0];
  yb1 = dboxSM2[1] - dboxUM2[1];
  yb2 = yb1 - dboxUM2[1] - 0.1 -  dboxUM2[1];
  yb3 = yb2 - dboxUM2[1] - 0.1 -  dboxUM2[1];


  //PRESHOWER  
  gMC->Gspos("EUM2", 1, "ESMB", xb1, yb1, 0., 0, "ONLY");
  gMC->Gspos("EUM2", 2, "ESMB", xb2, yb1, 0., 0, "ONLY");
  gMC->Gspos("EUM2", 3, "ESMB", xb1, yb2, 0., 0, "ONLY");
  gMC->Gspos("EUM2", 4, "ESMB", xb2, yb2, 0., 0, "ONLY");
  gMC->Gspos("EUM2", 5, "ESMB", xb1, yb3, 0., 0, "ONLY");
  gMC->Gspos("EUM2", 6, "ESMB", xb2, yb3, 0., 0, "ONLY");
  
  //VETO
  gMC->Gspos("EUV2", 1, "EMVB", xb1, yb1, 0., 0, "ONLY");
  gMC->Gspos("EUV2", 2, "EMVB", xb2, yb1, 0., 0, "ONLY");
  gMC->Gspos("EUV2", 3, "EMVB", xb1, yb2, 0., 0, "ONLY");
  gMC->Gspos("EUV2", 4, "EMVB", xb2, yb2, 0., 0, "ONLY");
  gMC->Gspos("EUV2", 5, "EMVB", xb1, yb3, 0., 0, "ONLY");
  gMC->Gspos("EUV2", 6, "EMVB", xb2, yb3, 0., 0, "ONLY");
  
  // Make smiliar stucture for lead as for PMD plane
  //================================================

  // 2 types of Rectangular shaped supermodules (BOX) 
  //each with 6 unit modules 
  
  // volume for SUPERMODULE ESMPbA 
  //Space added to provide a gapping for HV between UM's

  Float_t dboxSMPb1[3];
  dboxSMPb1[0] = 3.0*dboxUM1[0] + (2.0*0.075);
  dboxSMPb1[1] = 2.0*dboxUM1[1] +  0.05;
  dboxSMPb1[2] = fgkThLead/2.;
  
  gMC->Gsvolu("ESPA","BOX", idtmed[698], dboxSMPb1, 3);
  gMC->Gsatt("ESPA", "SEEN", 0);
  

  //Position the 6 unit modules in ESMPbA
  Float_t xpa1,xpa2,xpa3,ypa1,ypa2; 
  xpa1 = -dboxSMPb1[0] + dboxUM1[0];
  xpa2 = xpa1 + dboxUM1[0] + 0.15 + dboxUM1[0];
  xpa3 = xpa2 + dboxUM1[0] + 0.15 + dboxUM1[0];
  ypa1 = dboxSMPb1[1]  - dboxUM1[1];
  ypa2 = ypa1 - dboxUM1[1] - 0.1 - dboxUM1[1];


  gMC->Gspos("EPB1", 1, "ESPA", xpa1, ypa1, 0., 0, "ONLY");
  gMC->Gspos("EPB1", 2, "ESPA", xpa2, ypa1, 0., 0, "ONLY");
  gMC->Gspos("EPB1", 3, "ESPA", xpa3, ypa1, 0., 0, "ONLY");
  gMC->Gspos("EPB1", 4, "ESPA", xpa1, ypa2, 0., 0, "ONLY");
  gMC->Gspos("EPB1", 5, "ESPA", xpa2, ypa2, 0., 0, "ONLY");
  gMC->Gspos("EPB1", 6, "ESPA", xpa3, ypa2, 0., 0, "ONLY");


  // volume for SUPERMODULE ESMPbB 
  //Space is added to provide a gapping for HV between UM's
  Float_t dboxSMPb2[3];
  dboxSMPb2[0] = 2.0*dboxUM2[0] + 0.075;
  dboxSMPb2[1] = 3.0*dboxUM2[1] + (2.0*0.05);
  dboxSMPb2[2] = fgkThLead/2.;

  gMC->Gsvolu("ESPB","BOX", idtmed[698], dboxSMPb2, 3);
  gMC->Gsatt("ESPB", "SEEN", 0);
 
  //Position the 6 unit modules in ESMPbB
  Float_t xpb1,xpb2,ypb1,ypb2,ypb3; 
  xpb1 = -dboxSMPb2[0] + dboxUM2[0];
  xpb2 = xpb1 + dboxUM2[0] + 0.15 + dboxUM2[0];
  ypb1 = dboxSMPb2[1]  - dboxUM2[1];
  ypb2 = ypb1 - dboxUM2[1] - 0.1 -  dboxUM2[1];
  ypb3 = ypb2 - dboxUM2[1] - 0.1 -  dboxUM2[1];


  gMC->Gspos("EPB2", 1, "ESPB", xpb1, ypb1, 0., 0, "ONLY");
  gMC->Gspos("EPB2", 2, "ESPB", xpb2, ypb1, 0., 0, "ONLY");
  gMC->Gspos("EPB2", 3, "ESPB", xpb1, ypb2, 0., 0, "ONLY");
  gMC->Gspos("EPB2", 4, "ESPB", xpb2, ypb2, 0., 0, "ONLY");
  gMC->Gspos("EPB2", 5, "ESPB", xpb1, ypb3, 0., 0, "ONLY");
  gMC->Gspos("EPB2", 6, "ESPB", xpb2, ypb3, 0., 0, "ONLY");


  //---------------------------------------------------
  /// ALICE PMD FEE BOARDS IMPLEMENTATION
  // Dt: 25th February 2006 
  // - M.M.  Mondal, S.K. Prasad and P.K. Netrakanti
  //---------------------------------------------------

  //FEE boards 
  // It is FR4 board of length 7cm
  // breadth of 2.4 cm and thickness 0.1cm
  Float_t dboxFEE[3];
  dboxFEE[0] = 0.05;
  dboxFEE[1] = 3.50;
  dboxFEE[2] = 1.20;

  gMC->Gsvolu("EFEE","BOX", idtmed[607], dboxFEE, 3);
  gMC->Gsatt("EFEE", "SEEN", 0);
  gMC->Gsatt("EFEE", "COLO", 4);

  //Mother volume to accomodate FEE boards
  // It should have the dimension 
  // as the back plane or the 
  //corresponding UM
  //TYPE A
  //------------------------------------------------------//

  Float_t dboxFEEBPlaneA[3];
  dboxFEEBPlaneA[0]   = dboxBPlaneA[0]; //dbox3[0]; 
  dboxFEEBPlaneA[1]   = dboxBPlaneA[1];//dbox3[1];       
  dboxFEEBPlaneA[2]   = 1.2;
  
  //Volume of same dimension as Back PLane of Material AIR
  gMC->Gsvolu("EFBA","BOX", idtmed[698], dboxFEEBPlaneA, 3);
  gMC->Gsatt("EFBA", "SEEN", 0);

  //TYPE B
  Float_t dboxFEEBPlaneB[3];
  dboxFEEBPlaneB[0]   = dboxBPlaneB[0]; //dbox4[0]; 
  dboxFEEBPlaneB[1]   = dboxBPlaneB[1];//dbox4[1];       
  dboxFEEBPlaneB[2]   = 1.2;
  
  //Back PLane PCB of MAterial G10
  gMC->Gsvolu("EFBB","BOX", idtmed[698], dboxFEEBPlaneB, 3);
  gMC->Gsatt("EFBB", "SEEN", 0);

  //Placing the FEE boards in the Mother volume of AIR

  //Type A 

  Float_t xFee; // X-position of FEE board
  Float_t yFee; // Y-position of FEE board
  Float_t zFee = 0.0; // Z-position of FEE board

  Float_t xA    = 0.25; //distance from the border to 1st FEE board
  Float_t yA    = 4.00; //distance from the border to 1st FEE board
  Float_t xSepa = 1.70; //Distance between two FEE boards
  Float_t ySepa = 8.00; //Distance between two FEE boards

  
  // FEE Boards EFEE placed inside EFBA
  number = 1;
  yFee =  dboxFEEBPlaneA[1] - yA;  
  for (i = 1; i <= 6; ++i) 
    {
      xFee = -dboxFEEBPlaneA[0] + xA; 
      for (j = 1; j <= 12; ++j) 
	{
	  gMC->Gspos("EFEE", number, "EFBA", xFee,yFee,zFee, 0, "ONLY");
	  xFee += xSepa;
	  number += 1;
	}
      yFee -= ySepa; 
    }
  // FEE Boards EFEE placed inside EFBB
  number = 1;
  yFee =  dboxFEEBPlaneB[1] - yA;  
  for (i = 1; i <= 3; ++i) 
    {
      xFee = -dboxFEEBPlaneB[0] + xA; 
      for (j = 1; j <= 24; ++j) 
	{
	  gMC->Gspos("EFEE", number, "EFBB", xFee,yFee,zFee, 0, "ONLY");
	  xFee += xSepa;
	  number += 1;
	}
      yFee -= ySepa; 
    }


  //Distance between the two backplanes of two UMs
  //in x-direction is 0.92 and ydirection is 0.95
  Float_t dboxEFSA[3];
  dboxEFSA[0] = 3.0*dboxFEEBPlaneA[0] + 0.92;
  dboxEFSA[1] = 2.0*dboxFEEBPlaneA[1] + (0.95/2.0);
  dboxEFSA[2] = dboxFEEBPlaneA[2];

  //Type A
  gMC->Gsvolu("EFSA","BOX", idtmed[698],dboxEFSA, 3);
  gMC->Gsatt("EFSA", "SEEN", 0);

  //Distance between the two backplanes of two UMs
  //in x-direction is 0.92 and ydirection is 0.95
  Float_t dboxEFSB[3];
  dboxEFSB[0] = 2.0*dboxFEEBPlaneB[0] + (0.938/2.0);
  dboxEFSB[1] = 3.0*dboxFEEBPlaneB[1] + 1.05;
  dboxEFSB[2] = dboxFEEBPlaneB[2];

  //Type A
  gMC->Gsvolu("EFSB","BOX", idtmed[698],dboxEFSB, 3);
  gMC->Gsatt("EFSB", "SEEN", 0);


  Float_t xfs1,xfs2,xfs3,yfs1,yfs2,yfs3; 
  xfs1 = -dboxEFSA[0] + dboxFEEBPlaneA[0];
  xfs2 = xfs1 + dboxFEEBPlaneA[0] +  0.92 + dboxFEEBPlaneA[0];
  xfs3 = xfs2 + dboxFEEBPlaneA[0] +  0.92 + dboxFEEBPlaneA[0];
  yfs1 = dboxEFSA[1] - dboxFEEBPlaneA[1];
  yfs2 = yfs1 - dboxFEEBPlaneA[1] - 0.95 - dboxFEEBPlaneA[1];



  gMC->Gspos("EFBA", 1, "EFSA", xfs1, yfs1, 0., 0, "ONLY");
  gMC->Gspos("EFBA", 2, "EFSA", xfs2, yfs1, 0., 0, "ONLY");
  gMC->Gspos("EFBA", 3, "EFSA", xfs3, yfs1, 0., 0, "ONLY");
  gMC->Gspos("EFBA", 4, "EFSA", xfs1, yfs2, 0., 0, "ONLY");
  gMC->Gspos("EFBA", 5, "EFSA", xfs2, yfs2, 0., 0, "ONLY");
  gMC->Gspos("EFBA", 6, "EFSA", xfs3, yfs2, 0., 0, "ONLY");


  //Type B positioning

  xfs1 = -dboxEFSB[0] + dboxFEEBPlaneB[0];
  xfs2 = xfs1 + dboxFEEBPlaneB[0] + 0.938 + dboxFEEBPlaneB[0];
  yfs1 = dboxEFSB[1] - dboxFEEBPlaneB[1];
  yfs2 = yfs1 - dboxFEEBPlaneB[1] - 1.05 - dboxFEEBPlaneB[1];
  yfs3 = yfs2 - dboxFEEBPlaneB[1] - 1.05 - dboxFEEBPlaneB[1];



  gMC->Gspos("EFBB", 1, "EFSB", xfs1, yfs1, 0., 0, "ONLY");
  gMC->Gspos("EFBB", 2, "EFSB", xfs2, yfs1, 0., 0, "ONLY");
  gMC->Gspos("EFBB", 3, "EFSB", xfs1, yfs2, 0., 0, "ONLY");
  gMC->Gspos("EFBB", 4, "EFSB", xfs2, yfs2, 0., 0, "ONLY");
  gMC->Gspos("EFBB", 5, "EFSB", xfs1, yfs3, 0., 0, "ONLY");
  gMC->Gspos("EFBB", 6, "EFSB", xfs2, yfs3, 0., 0, "ONLY");


}
 
//_____________________________________________________________________________

void AliPMDv1::CreatePMD()
{
  //
  // Create final detector from supermodules
  // -- Author : Bedanga and Viyogi June 2003

  Float_t   zp;
  Int_t jhrot12,jhrot13, irotdm;
  Int_t *idtmed = fIdtmed->GetArray()-599;
  
  //VOLUMES Names : begining with "E" for all PMD volumes, 

  // --- DEFINE Iron volumes  for SM A
  //   Fe Support 
  Float_t dboxFea[3];
  dboxFea[0] = fSMLengthax;
  dboxFea[1] = fSMLengthay;
  dboxFea[2] = fgkThSteel/2.;
  
  gMC->Gsvolu("EFEA","BOX", idtmed[618], dboxFea, 3);
  gMC->Gsatt ("EFEA", "SEEN", 0);

  // --- DEFINE Iron volumes  for SM B
  
  //   Fe Support 
  Float_t dboxFeb[3];
  dboxFeb[0] = fSMLengthbx;
  dboxFeb[1] = fSMLengthby;
  dboxFeb[2] = fgkThSteel/2.;
  
  gMC->Gsvolu("EFEB","BOX", idtmed[618], dboxFeb, 3);
  gMC->Gsatt ("EFEB", "SEEN", 0);

  AliMatrix(irotdm, 90., 0.,  90.,  90., 180., 0.);
  AliMatrix(jhrot12, 90., 180., 90., 270., 0., 0.);
  AliMatrix(jhrot13, 90., 240., 90., 330., 0., 0.);

  // Gaspmd, the dimension of RECTANGULAR mother volume of PMD,
  // Four mother volumes EPM1,EPM2 for A-type and 
  // volumes EPM3 and EPM4 for B-type. Four to create a hole
  // and avoid overlap with beam pipe

  Float_t gaspmd[3];
  gaspmd[0] = fSMLengthax;
  gaspmd[1] = fSMLengthay;
  gaspmd[2] = fSMthick;

  gMC->Gsvolu("EPM1", "BOX", idtmed[698], gaspmd, 3);
  gMC->Gsatt("EPM1", "SEEN", 1);
  gMC->Gsvolu("EPM2", "BOX", idtmed[698], gaspmd, 3);
  gMC->Gsatt("EPM2", "SEEN", 1);

  //Complete detector for Type A
  //Position Super modules type A for both CPV and PMD in EPMD  
  Float_t zpsa,zpba,zfea,zcva,zfee; 

  // zpsa = - gaspmd[2] + fSMthick/2.;
  // -2.5 is given to place PMD at -361.5 
  // BM : In future after putting proper electronics
  // -2.5 will be replaced by -gaspmd[2]

  //TYPE A
  //Fee board
  zfee=-gaspmd[2] + 1.2;
  gMC->Gspos("EFSA", 1, "EPM1", 0., 0., zfee, 0, "ONLY");
  gMC->Gspos("EFSA", 2, "EPM2", 0., 0., zfee, jhrot12, "ONLY");
  //VETO
  zcva = zfee + 1.2 + fDthick;
  gMC->Gspos("EMVA", 1, "EPM1", 0., 0., zcva, 0, "ONLY");
  gMC->Gspos("EMVA", 2, "EPM2", 0., 0., zcva, jhrot12, "ONLY");
  //Iron support
  zfea = zcva + fDthick + fgkThSteel/2.;
  gMC->Gspos("EFEA", 1, "EPM1", 0., 0., zfea, 0, "ONLY");
  gMC->Gspos("EFEA", 2, "EPM2", 0., 0., zfea, 0, "ONLY");
  //Lead
  zpba=zfea+fgkThSteel/2.+ fgkThLead/2.;
  gMC->Gspos("ESPA", 1, "EPM1", 0., 0., zpba, 0, "ONLY");
  gMC->Gspos("ESPA", 2, "EPM2", 0., 0., zpba, 0, "ONLY");
  //Preshower
  zpsa = zpba + fgkThLead/2. + fDthick;
  gMC->Gspos("ESMA", 1, "EPM1", 0., 0., zpsa, 0, "ONLY");
  gMC->Gspos("ESMA", 2, "EPM2", 0., 0., zpsa, jhrot12, "ONLY");
  //FEE boards
  zfee=zpsa + fDthick + 1.2;
  gMC->Gspos("EFSA", 3, "EPM1", 0., 0., zfee, 0, "ONLY");
  gMC->Gspos("EFSA", 4, "EPM2", 0., 0., zfee, jhrot12, "ONLY");

 
  //TYPE - B
  gaspmd[0] = fSMLengthbx; 
  gaspmd[1] = fSMLengthby; 
  gaspmd[2] = fSMthick; 

  gMC->Gsvolu("EPM3", "BOX", idtmed[698], gaspmd, 3);
  gMC->Gsatt("EPM3", "SEEN", 1);
  gMC->Gsvolu("EPM4", "BOX", idtmed[698], gaspmd, 3);
  gMC->Gsatt("EPM4", "SEEN", 1);

  //Complete detector for Type B
  //Position Super modules type B for both CPV and PMD in EPMD  
  Float_t zpsb,zpbb,zfeb,zcvb; 
  // zpsb = - gaspmd[2] + fSMthick/2.;
  // -2.5 is given to place PMD at -361.5 
  // BM: In future after putting proper electronics
  // -2.5 will be replaced by -gaspmd[2]

 //Fee board
  zfee=-gaspmd[2] + 1.2;
  gMC->Gspos("EFSB", 5, "EPM3", 0., 0., zfee, 0, "ONLY");
  gMC->Gspos("EFSB", 6, "EPM4", 0., 0., zfee, jhrot12, "ONLY");
  //VETO
  zcvb= zfee + 1.2 + fDthick;
  gMC->Gspos("EMVB", 3, "EPM3", 0., 0., zcvb, 0, "ONLY");
  gMC->Gspos("EMVB", 4, "EPM4", 0., 0., zcvb, jhrot12, "ONLY");

  //IRON SUPPORT
  zfeb= zcvb + fDthick +  fgkThSteel/2.;
  gMC->Gspos("EFEB", 3, "EPM3", 0., 0., zfeb, 0, "ONLY");
  gMC->Gspos("EFEB", 4, "EPM4", 0., 0., zfeb, 0, "ONLY");
  //LEAD
  zpbb= zfeb + fgkThSteel/2.+ fgkThLead/2.;
  gMC->Gspos("ESPB", 3, "EPM3", 0., 0., zpbb, 0, "ONLY");
  gMC->Gspos("ESPB", 4, "EPM4", 0., 0., zpbb, 0, "ONLY");
  //PRESHOWER
  zpsb = zpbb + fgkThLead/2.+ fDthick;
  gMC->Gspos("ESMB", 3, "EPM3", 0., 0., zpsb, 0, "ONLY");
  gMC->Gspos("ESMB", 4, "EPM4", 0., 0., zpsb, jhrot12, "ONLY");
  //FEE boards
  zfee=zpsb + fDthick + 1.2;
  gMC->Gspos("EFSB", 7, "EPM3", 0., 0., zfee, 0, "ONLY");
  gMC->Gspos("EFSB", 8, "EPM4", 0., 0., zfee, jhrot12, "ONLY");


  // --- Place the EPMD in ALICE 
  //Z-distance of PMD from Interaction Point
  zp = fgkZdist;

  //X and Y-positions of the PMD planes
  Float_t xfinal,yfinal; 
  Float_t xsmb,ysmb;
  Float_t xsma,ysma;

  xfinal = fSMLengthax + 0.48/2 + fSMLengthbx;
  yfinal = fSMLengthay + 0.20/2 + fSMLengthby;
  

  xsma =  xfinal  - fSMLengthax;
  ysma =  yfinal  - fSMLengthay;
  xsmb =  -xfinal + fSMLengthbx;
  ysmb =  yfinal  - fSMLengthby;


//Position Full PMD in ALICE   
//
//   EPM1      EPM3
//
//   EPM4      EPM2
// (rotated   (rotated EPM1)
//  EPM3)
//
  gMC->Gspos("EPM1", 1, "ALIC",  xsma,ysma,zp,  0, "ONLY");
  gMC->Gspos("EPM2", 1, "ALIC", -xsma,-ysma,zp, 0, "ONLY");
  gMC->Gspos("EPM3", 1, "ALIC",  xsmb,ysmb,zp,  0, "ONLY");
  gMC->Gspos("EPM4", 1, "ALIC", -xsmb,-ysmb,zp, 0, "ONLY");
}

 
//_____________________________________________________________________________
void AliPMDv1::DrawModule() const
{
  // Draw a shaded view of the Photon Multiplicity Detector
  //
  //  cout << " Inside Draw Modules " << endl;

  gMC->Gsatt("*", "seen", -1);
  gMC->Gsatt("alic", "seen", 0);
  //
  // Set the visibility of the components
  // 
  gMC->Gsatt("ECAR","seen",0);
  gMC->Gsatt("ECCU","seen",1);
  gMC->Gsatt("EST1","seen",1);
  gMC->Gsatt("EST2","seen",1);
  gMC->Gsatt("EUM1","seen",1);
  gMC->Gsatt("EUM2","seen",1);
  gMC->Gsatt("ESMA","seen",1);
  gMC->Gsatt("EPMD","seen",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 3000, -3000, 3000, -6000, 6000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 22, 20.5, .02, .02);
  gMC->Gdhead(1111, "Photon Multiplicity Detector Version 1");

  //gMC->Gdman(17, 5, "MAN");
  gMC->Gdopt("hide", "off");

  AliDebug(1,"Outside Draw Modules");
}

//_____________________________________________________________________________
void AliPMDv1::CreateMaterials()
{
  // Create materials for the PMD
  //
  // ORIGIN    : Y. P. VIYOGI 
  //
  //  cout << " Inside create materials " << endl;

  Int_t *idtmed = fIdtmed->GetArray()-599;
  Int_t isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  
  // --- Define the various materials for GEANT --- 

  AliMaterial(1, "Pb    $", 207.19, 82., 11.35, .56, 18.5);
  
  // Argon

  Float_t dAr   = 0.001782;   // --- Ar density in g/cm3 --- 
  Float_t x0Ar = 19.55 / dAr;
  AliMaterial(2, "Argon$", 39.95, 18., dAr, x0Ar, 6.5e4);

  // --- CO2 --- 

  Float_t aCO2[2] = { 12.,16. };
  Float_t zCO2[2] = { 6.,8. };
  Float_t wCO2[2] = { 1.,2. };
  Float_t dCO2    = 0.001977;
  AliMixture(3, "CO2  $", aCO2, zCO2, dCO2, -2, wCO2);

  AliMaterial(4, "Al   $", 26.98, 13., 2.7, 8.9, 18.5);

  // ArCO2

  Float_t aArCO2[3] = {39.948,12.0107,15.9994};
  Float_t zArCO2[3] = {18.,6.,8.};
  Float_t wArCO2[3] = {0.7,0.08,0.22};
  Float_t dArCO2    = dAr * 0.7 + dCO2 * 0.3;
  AliMixture(5, "ArCO2$", aArCO2, zArCO2, dArCO2, 3, wArCO2);

  AliMaterial(6, "Fe   $", 55.85, 26., 7.87, 1.76, 18.5);

  // G10
  
  Float_t aG10[4]={1.,12.011,15.9994,28.086};
  Float_t zG10[4]={1.,6.,8.,14.};
  Float_t wG10[4]={0.15201,0.10641,0.49444,0.24714};
  AliMixture(8,"G10",aG10,zG10,1.7,4,wG10);
  
  AliMaterial(15, "Cu   $", 63.54, 29., 8.96, 1.43, 15.);

  // Steel
  Float_t aSteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zSteel[4] = { 26.,24.,28.,14. };
  Float_t wSteel[4] = { .715,.18,.1,.005 };
  Float_t dSteel    = 7.88;
  AliMixture(19, "STAINLESS STEEL$", aSteel, zSteel, dSteel, 4, wSteel); 

  //Air

  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir1 = 1.20479E-10;
  Float_t dAir = 1.20479E-3;
  AliMixture(98, "Vacum$", aAir,  zAir, dAir1, 4, wAir);
  AliMixture(99, "Air  $", aAir,  zAir, dAir , 4, wAir);

  // Define tracking media 
  AliMedium(1,  "Pb conv.$", 1,  0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
  AliMedium(4,  "Al      $", 4,  0, 0, isxfld, sxmgmx, .1, .1, .01, .1);
  AliMedium(5,  "ArCO2   $", 5,  1, 0, isxfld, sxmgmx, .1, .1, .10, .1);
  AliMedium(6,  "Fe      $", 6,  0, 0, isxfld, sxmgmx, .1, .1, .01, .1);
  AliMedium(8,  "G10plate$", 8,  0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
  AliMedium(15, "Cu      $", 15, 0, 0, isxfld, sxmgmx, .1, .1, .01, .1);
  AliMedium(19, "S  steel$", 19, 0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
  AliMedium(98, "Vacuum  $", 98, 0, 0, isxfld, sxmgmx, 1., .1, .10, 10);
  AliMedium(99, "Air gaps$", 99, 0, 0, isxfld, sxmgmx, 1., .1, .10, .1);
  
  // --- Generate explicitly delta rays in the iron, aluminium and lead --- 
  gMC->Gstpar(idtmed[600], "LOSS", 3.);
  gMC->Gstpar(idtmed[600], "DRAY", 1.);
  
  gMC->Gstpar(idtmed[603], "LOSS", 3.);
  gMC->Gstpar(idtmed[603], "DRAY", 1.);
  
  gMC->Gstpar(idtmed[604], "LOSS", 3.);
  gMC->Gstpar(idtmed[604], "DRAY", 1.);
  
  gMC->Gstpar(idtmed[605], "LOSS", 3.);
  gMC->Gstpar(idtmed[605], "DRAY", 1.);
  
  gMC->Gstpar(idtmed[607], "LOSS", 3.);
  gMC->Gstpar(idtmed[607], "DRAY", 1.);
  
  // --- Energy cut-offs in the Pb and Al to gain time in tracking --- 
  // --- without affecting the hit patterns --- 
  gMC->Gstpar(idtmed[600], "CUTGAM", 1e-4);
  gMC->Gstpar(idtmed[600], "CUTELE", 1e-4);
  gMC->Gstpar(idtmed[600], "CUTNEU", 1e-4);
  gMC->Gstpar(idtmed[600], "CUTHAD", 1e-4);

  gMC->Gstpar(idtmed[605], "CUTGAM", 1e-4);
  gMC->Gstpar(idtmed[605], "CUTELE", 1e-4);
  gMC->Gstpar(idtmed[605], "CUTNEU", 1e-4);
  gMC->Gstpar(idtmed[605], "CUTHAD", 1e-4);

  gMC->Gstpar(idtmed[603], "CUTGAM", 1e-4);
  gMC->Gstpar(idtmed[603], "CUTELE", 1e-4);
  gMC->Gstpar(idtmed[603], "CUTNEU", 1e-4);
  gMC->Gstpar(idtmed[603], "CUTHAD", 1e-4);
//   gMC->Gstpar(idtmed[609], "CUTGAM", 1e-4);
//   gMC->Gstpar(idtmed[609], "CUTELE", 1e-4);
//   gMC->Gstpar(idtmed[609], "CUTNEU", 1e-4);
//   gMC->Gstpar(idtmed[609], "CUTHAD", 1e-4);
  // --- Prevent particles stopping in the gas due to energy cut-off --- 
  gMC->Gstpar(idtmed[604], "CUTGAM", 1e-5);
  gMC->Gstpar(idtmed[604], "CUTELE", 1e-5);
  gMC->Gstpar(idtmed[604], "CUTNEU", 1e-5);
  gMC->Gstpar(idtmed[604], "CUTHAD", 1e-5);
  gMC->Gstpar(idtmed[604], "CUTMUO", 1e-5);

  AliDebug(1,"Outside create materials");

}

//_____________________________________________________________________________
void AliPMDv1::Init()
{
  //
  // Initialises PMD detector after it has been built
  //

  //
  AliDebug(2,"Inside Init");
  AliDebug(2,"PMD simulation package (v1) initialised");
  AliDebug(2,"parameters of pmd");
  AliDebug(2,Form("%10.2f %10.2f %10.2f %10.2f\n",
		  fgkCellRadius,fgkCellWall,fgkCellDepth,fgkZdist));
  Int_t *idtmed = fIdtmed->GetArray()-599;
  fMedSens=idtmed[605-1];

}

//_____________________________________________________________________________
void AliPMDv1::StepManager()
{
  //
  // Called at each step in the PMD
  //

  Int_t   copy;
  Float_t hits[4], destep;
  Float_t center[3] = {0,0,0};
  Int_t   vol[10];
  //  const char *namep;
  
  if(gMC->CurrentMedium() == fMedSens && (destep = gMC->Edep())) {
  
    gMC->CurrentVolID(copy);
    //     namep=gMC->CurrentVolName();
    // printf("Current vol  is %s \n",namep);
    vol[0]=copy;

    gMC->CurrentVolOffID(1,copy);
    //namep=gMC->CurrentVolOffName(1);
    // printf("Current vol 11 is %s \n",namep);
    vol[1]=copy;

    gMC->CurrentVolOffID(2,copy);
    //namep=gMC->CurrentVolOffName(2);
    //printf("Current vol 22 is %s \n",namep);
    vol[2]=copy;

    //	if(strncmp(namep,"EHC1",4))vol[2]=1;

    gMC->CurrentVolOffID(3,copy);
    // namep=gMC->CurrentVolOffName(3);
    //printf("Current vol 33 is %s \n",namep);
    vol[3]=copy;

    gMC->CurrentVolOffID(4,copy);
    // namep=gMC->CurrentVolOffName(4);
    // printf("Current vol 44 is %s \n",namep);
    vol[4]=copy;

    gMC->CurrentVolOffID(5,copy);
    // namep=gMC->CurrentVolOffName(5);
    // printf("Current vol 55 is %s \n",namep);
    vol[5]=copy;

    gMC->CurrentVolOffID(6,copy);
    // namep=gMC->CurrentVolOffName(6);
    // printf("Current vol 66 is %s \n",namep);
    vol[6]=copy;

    gMC->CurrentVolOffID(7,copy);
    //  namep=gMC->CurrentVolOffName(7);
    // printf("Current vol 77 is %s \n",namep);
    vol[7]=copy;

    gMC->CurrentVolOffID(8,copy);
    // namep=gMC->CurrentVolOffName(8);
    // printf("Current vol 88 is %s \n",namep);
    vol[8]=copy;


    gMC->CurrentVolOffID(9,copy);
    // namep=gMC->CurrentVolOffName(9);
    // printf("Current vol 99 is %s \n",namep);
    vol[9]=copy;


    // printf("volume number %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %10.3f \n",vol[0],vol[1],vol[2],vol[3],vol[4],vol[5],vol[6],vol[7],vol[8],vol[9],destep*1000000);
    
    gMC->Gdtom(center,hits,1);
    hits[3] = destep*1e9; //Number in eV
    AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);

  }
}

  
//------------------------------------------------------------------------
// Get parameters

void AliPMDv1::GetParameters()
{
  // This gives all the parameters of the detector
  // such as Length of Supermodules, type A, type B,
  // thickness of the Supermodule
  //
  
  fSMLengthax = 32.7434;
  //The total length in X is due to the following components
  // Factor 3 is because of 3 module length in X for this type
  // fgkNcolUM1*fgkCellRadius (48 x 0.25): Total span of each module in X
  // fgkCellRadius/2. : There is offset of 1/2 cell
  // 0.05+0.05 : Insulation gaps etc
  // fgkSSBoundary (0.3) : Boundary frame
  // double XA = 3.0*((fgkCellRadius/fgkSqroot3by2*fgkNcolUM1)-(fgkCellRadius*fgkSqroot3*(fgkNcolUM1-1)/6.)+(2.0*fgkGap)+(2.0*fgkGap)+fgkSSBoundary) + (2.0*0.075);

  fSMLengthbx = 42.5886;
  //The total length in X is due to the following components
  // Factor 2 is because of 2 module length in X for this type
  // fgkNcolUM2*fgkCellRadius (96 x 0.25): Total span of each module in X
  // fgkCellRadius/2. : There is offset of 1/2 cell
  // 0.05+0.05 : Insulation gaps etc
  // fgkSSBoundary (0.3) : Boundary frame
  //double XB = 2.0*((fgkCellRadius/fgkSqroot3by2*fgkNcolUM2)-(fgkCellRadius*fgkSqroot3*(fgkNcolUM2-1)/6.)+(2.0*fgkGap)+(2.0*fgkGap)+fgkSSBoundary) + 0.075; 



  fSMLengthay = 49.1;
  //The total length in Y is due to the following components
  // Factor 2 is because of 2 module length in Y for this type
  // fgkCellRadius/fgkSqroot3by2)*fgkNrowUM1 (0.25/sqrt3/2 * 96): Total span of each module in Y
  //  of strips
  // 0.05+0.05 : Insulation gaps etc
  // fgkSSBoundary (0.3) : Boundary frame
  // double  YA = 2.0*(fgkNrowUM1*fgkCellRadius+fgkCellRadius/2.+(2.0*fgkGap)+(2.0*fgkGap)+fgkSSBoundary) +  0.05;

  fSMLengthby =  37.675;
  //The total length in Y is due to the following components
  // Factor 3 is because of 3 module length in Y for this type
  // fgkCellRadius/fgkSqroot3by2)*fgkNrowUM2 (0.25/sqrt3/2 * 48): Total span of each module in Y
  //  of strips
  // 0.05+0.05 : Insulation gaps etc
  // fgkSSBoundary (0.3) : Boundary frame
    //double YB = 3.0*((fgkNrowUM2*fgkCellRadius + fgkCellRadius/2.)+(2.0*fgkGap)+(2.0*fgkGap)+fgkSSBoundary) + (2.0*0.05);


  //Thickness of a pre/veto plane 
  fDthick     = fgkThSS/2. +0.15;

  //Thickness of the PMD ; 2.4 added for FEE boards 
    fSMthick    = 2.0*(fgkThSS/2. +0.15)
                +fgkThSteel/2.+fgkThLead/2.0 + 2.4;


  
}
// ---------------------------------------------------------------
void AliPMDv1::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  // 
  SetSectorAlignable();

}
// ----------------------------------------------------------------
void AliPMDv1::SetSectorAlignable() const
{
  //

  TString vpsector = "ALIC_1/EPM";
  TString vpappend = "_1";

  TString snsector="PMD/Sector";

  TString volpath, symname;
  
  for(Int_t cnt=1; cnt<=4; cnt++){
    volpath = vpsector;
    volpath += cnt;
    volpath += vpappend;
    symname = snsector;
    symname += cnt;
    gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data());
  }
}
// ------------------------------------------------------------------
