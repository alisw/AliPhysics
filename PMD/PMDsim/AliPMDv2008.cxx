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
/* $Id: AliPMDv1.cxx 18594 2007-05-15 13:28:06Z hristov $ */

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

#include <Riostream.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>
#include <TVirtualMC.h>

#include "AliConst.h" 
#include "AliLog.h"
#include "AliMC.h"
#include "AliMagF.h" 
#include "AliPMDv2008.h"
#include "AliRun.h"

const Int_t   AliPMDv2008::fgkNcolUM1    = 48;  // Number of cols in UM, type 1
const Int_t   AliPMDv2008::fgkNcolUM2    = 96;  // Number of cols in UM, type 2
const Int_t   AliPMDv2008::fgkNrowUM1    = 96;  // Number of rows in UM, type 1
const Int_t   AliPMDv2008::fgkNrowUM2    = 48;  // Number of rows in UM, type 2
const Float_t AliPMDv2008::fgkCellRadius = 0.25;     // Radius of a hexagonal cell
const Float_t AliPMDv2008::fgkCellWall   = 0.02;     // Thickness of cell Wall
const Float_t AliPMDv2008::fgkCellDepth  = 0.50;     // Gas thickness
const Float_t AliPMDv2008::fgkThBase     = 0.2;      // Thickness of Base plate
const Float_t AliPMDv2008::fgkThBKP      = 0.1;      // Thickness of Back plane
const Float_t AliPMDv2008::fgkThAir      = 1.03;      // Thickness of Air
const Float_t AliPMDv2008::fgkThPCB      = 0.16;     // Thickness of PCB
const Float_t AliPMDv2008::fgkThLead     = 1.5;      // Thickness of Pb
const Float_t AliPMDv2008::fgkThSteel    = 0.5;      // Thickness of Steel
const Float_t AliPMDv2008::fgkGap        = 0.025;    // Air Gap
const Float_t AliPMDv2008::fgkZdist      = 361.5;    // z-position of the detector
const Float_t AliPMDv2008::fgkSqroot3    = 1.7320508;// Square Root of 3
const Float_t AliPMDv2008::fgkSqroot3by2 = 0.8660254;// Square Root of 3 by 2
const Float_t AliPMDv2008::fgkSSBoundary = 0.3;
const Float_t AliPMDv2008::fgkThSS       = 1.03;
const Float_t AliPMDv2008::fgkThG10      = 1.03;
ClassImp(AliPMDv2008)
 
//_____________________________________________________________________________
AliPMDv2008::AliPMDv2008():
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
AliPMDv2008::AliPMDv2008(const char *name, const char *title):
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
void AliPMDv2008::CreateGeometry()
{
  // Create geometry for Photon Multiplicity Detector

  GetParameters();
  CreateSupermodule();
  CreatePMD();
}

//_____________________________________________________________________________
void AliPMDv2008::CreateSupermodule()
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
  
  TVirtualMC::GetMC()->Gsvolu("ECAR", "PGON", idtmed[604], hexd2,10);
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

  TVirtualMC::GetMC()->Gsvolu("ECCU", "PGON", idtmed[614], hexd1,10);

  // Place  inner hex (sensitive volume) inside outer hex (copper)
  
  TVirtualMC::GetMC()->Gspos("ECAR", 1, "ECCU", 0., 0., 0., 0, "ONLY");
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
  
  TVirtualMC::GetMC()->Gsvolu("EST1","BOX", idtmed[698], dbox1, 3);

  // volume for second strip EST2 


  Float_t dbox2[3];
  dbox2[1] = fgkNrowUM2*fgkCellRadius;
  dbox2[0] = dbox1[0];
  dbox2[2] = dbox1[2];

  TVirtualMC::GetMC()->Gsvolu("EST2","BOX", idtmed[698], dbox2, 3);

  // Place hexagonal cells ECCU placed inside EST1 
  xb = 0.; 
  zb = 0.;
  yb = (dbox1[1]) - fgkCellRadius; 
  for (i = 1; i <= fgkNrowUM1; ++i) 
    {
      number = i;
      TVirtualMC::GetMC()->Gspos("ECCU", number, "EST1", xb,yb,zb, 0, "ONLY");
      yb -= (fgkCellRadius*2.);
    }

  // Place hexagonal cells ECCU placed inside EST2 
  xb = 0.; 
  zb = 0.;
  yb = (dbox2[1]) - fgkCellRadius; 
  for (i = 1; i <= fgkNrowUM2; ++i) 
    {
      number = i;
      TVirtualMC::GetMC()->Gspos("ECCU", number, "EST2", xb,yb,zb, 0, "ONLY");
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
  TVirtualMC::GetMC()->Gsvolu("EHC1","BOX", idtmed[698], dbox3, 3);
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
      TVirtualMC::GetMC()->Gspos("EST1",number, "EHC1", xb, yb , 0. , 0, "MANY");
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
  TVirtualMC::GetMC()->Gsvolu("EHC2","BOX", idtmed[698], dbox4, 3);

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
    TVirtualMC::GetMC()->Gspos("EST2",number, "EHC2", xb, yb , 0. ,0, "MANY");
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
  TVirtualMC::GetMC()->Gsvolu("EPCA","BOX", idtmed[607], dboxPcbA, 3);
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
  TVirtualMC::GetMC()->Gsvolu("EBKA","BOX", idtmed[607], dboxBPlaneA, 3);
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
  TVirtualMC::GetMC()->Gsvolu("ECGA","BOX", idtmed[698], dboxAir3A, 3);

  //FOR VETO
  //Air gap is a BOX of Material Air
  TVirtualMC::GetMC()->Gsvolu("ECVA","BOX", idtmed[698], dboxAir3A, 3);
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
  TVirtualMC::GetMC()->Gsvolu("EDGA","BOX", idtmed[607], dboxGGA, 3);

  //FOR VETO
  //G10 BOX 
  TVirtualMC::GetMC()->Gsvolu("EDVA","BOX", idtmed[607], dboxGGA, 3);

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
  TVirtualMC::GetMC()->Gsvolu("ESSA","BOX", idtmed[618], dboxSS1, 3);

  //FOR VETO
  //Stainless Steel boundary - Material Stainless Steel
  TVirtualMC::GetMC()->Gsvolu("ESVA","BOX", idtmed[618], dboxSS1, 3);

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
  TVirtualMC::GetMC()->Gspos("EDVA", 1, "ESVA", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECVA", 1, "EDVA", 0., 0., 0., 0, "ONLY");

//FOR PRESHOWER
  TVirtualMC::GetMC()->Gspos("EDGA", 1, "ESSA", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECGA", 1, "EDGA", 0., 0., 0., 0, "ONLY");

 // now other components, using Bedanga's code, but changing the values.
  //Positioning Bottom PCB, Honey Comb abd Top PCB in AIR

  //For veto plane
  //Positioning the Bottom 0.16 cm PCB
  Float_t zbpcb = -dboxAir3A[2] + (2.0*fgkGap) + fgkThPCB/2.;
  TVirtualMC::GetMC()->Gspos("EPCA", 1, "ECVA", 0., 0., zbpcb, 0, "ONLY");
  //Positioning the Honey Comb 0.5 cm
  Float_t zhc = zbpcb + fgkThPCB/2. + fgkCellDepth/2.;
  TVirtualMC::GetMC()->Gspos("EHC1", 1, "ECVA", 0., 0., zhc, 0, "ONLY");
  //Positioning the Top PCB 0.16 cm
  Float_t ztpcb = zhc + fgkCellDepth/2 + fgkThPCB/2.;
  TVirtualMC::GetMC()->Gspos("EPCA", 2, "ECVA", 0., 0., ztpcb, 0, "ONLY");


  //For Preshower plane the ordering is reversed
  //Positioning the Bottom 0.16 cm PCB
  zbpcb = -dboxAir3A[2] + fgkThPCB + fgkThPCB/2.;
  TVirtualMC::GetMC()->Gspos("EPCA", 1, "ECGA", 0., 0., zbpcb, 0, "ONLY");
  //Positioning the Honey Comb 0.5 cm
  zhc = zbpcb + fgkThPCB/2. + fgkCellDepth/2.;
  TVirtualMC::GetMC()->Gspos("EHC1", 1, "ECGA", 0., 0., zhc, 0, "ONLY");
  //Positioning the Top PCB 0.16 cm
  ztpcb = zhc + fgkCellDepth/2 + fgkThPCB/2.;
  TVirtualMC::GetMC()->Gspos("EPCA", 2, "ECGA", 0., 0., ztpcb, 0, "ONLY");




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
  TVirtualMC::GetMC()->Gsvolu("EUM1","BOX", idtmed[698], dboxUM1, 3);
  //FOR VETO
  TVirtualMC::GetMC()->Gsvolu("EUV1","BOX", idtmed[698], dboxUM1, 3);

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
  TVirtualMC::GetMC()->Gsvolu("EBPA","BOX", idtmed[607], dboxBaseA, 3);
  //----------------------------------------------------//  

  //FOR VETO
  //- Placing of all components of UM in AIR BOX EUM1--//
  //(1)   FIRST PUT THE BASE PLATE
  Float_t zbaseplate = -dboxUM1[2] + fgkThBase/2.;
  TVirtualMC::GetMC()->Gspos("EBPA", 1, "EUV1", 0., 0., zbaseplate, 0, "ONLY");

  //(2)   NEXT PLACING the SS BOX 
  Float_t zss = zbaseplate + fgkThBase/2. + fgkThSS/2.;
  TVirtualMC::GetMC()->Gspos("ESVA", 1, "EUV1", 0., 0., zss, 0, "ONLY");
  
  // (3) Positioning the Backplane PCB 0.1 cm
  Float_t zbkp = zss + fgkThSS/2. + fgkThBKP/2.;
  TVirtualMC::GetMC()->Gspos("EBKA", 1, "EUV1", 0., 0., zbkp, 0, "ONLY");

  //FOR PRESHOWER
  // (3) Positioning the Backplane PCB 0.1 cm
  zbkp = -dboxUM1[2] + fgkThBKP/2.;
  TVirtualMC::GetMC()->Gspos("EBKA", 1, "EUM1", 0., 0., zbkp, 0, "ONLY");

  //(2)   NEXT PLACING the SS BOX 
  zss = zbkp + fgkThBKP/2. + fgkThSS/2.;
  TVirtualMC::GetMC()->Gspos("ESSA", 1, "EUM1", 0., 0., zss, 0, "ONLY");
  
  //(1)   FIRST PUT THE BASE PLATE
  zbaseplate = zss + fgkThSS/2 + fgkThBase/2.;
  TVirtualMC::GetMC()->Gspos("EBPA", 1, "EUM1", 0., 0., zbaseplate, 0, "ONLY");
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
  TVirtualMC::GetMC()->Gsvolu("EPCB","BOX", idtmed[607], dboxPcbB, 3);
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
  TVirtualMC::GetMC()->Gsvolu("EBKB","BOX", idtmed[607], dboxBPlaneB, 3);
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
  TVirtualMC::GetMC()->Gsvolu("ECGB","BOX", idtmed[698], dboxAir3B, 3);
  //VETO
  TVirtualMC::GetMC()->Gsvolu("ECVB","BOX", idtmed[698], dboxAir3B, 3);

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
  TVirtualMC::GetMC()->Gsvolu("EDGB","BOX", idtmed[607], dboxGGB, 3);
  //VETO
  TVirtualMC::GetMC()->Gsvolu("EDVB","BOX", idtmed[607], dboxGGB, 3);
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
  TVirtualMC::GetMC()->Gsvolu("ESSB","BOX", idtmed[618], dboxSS2, 3);
  //VETO
  TVirtualMC::GetMC()->Gsvolu("ESVB","BOX", idtmed[618], dboxSS2, 3);
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
  TVirtualMC::GetMC()->Gspos("EDGB", 1, "ESSB", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECGB", 1, "EDGB", 0., 0., 0., 0, "ONLY");
  //VETO
  TVirtualMC::GetMC()->Gspos("EDVB", 1, "ESVB", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECVB", 1, "EDVB", 0., 0., 0., 0, "ONLY");

 // now other components, using Bedang's code, but changing the values.
  //Positioning Bottom PCB, Honey Comb abd Top PCB in AIR

  //VETO
  //Positioning the Bottom 0.16 cm PCB
  Float_t zbpcb2 = -dboxAir3B[2] + (2.0*fgkGap) + fgkThPCB/2.;
  TVirtualMC::GetMC()->Gspos("EPCB", 1, "ECVB", 0., 0., zbpcb2, 0, "ONLY");
  //Positioning the Honey Comb 0.5 cm
  Float_t zhc2 = zbpcb2 + fgkThPCB/2. + fgkCellDepth/2.;
  TVirtualMC::GetMC()->Gspos("EHC2", 1, "ECVB", 0., 0., zhc2, 0, "ONLY");
  //Positioning the Top PCB 0.16 cm
  Float_t ztpcb2 = zhc2 + fgkCellDepth/2 + fgkThPCB/2.;
  TVirtualMC::GetMC()->Gspos("EPCB", 2, "ECVB", 0., 0., ztpcb2, 0, "ONLY");

  //PRESHOWER
  //For preshower plane the ordering is reversed
  //Positioning the Bottom 0.16 cm PCB
  zbpcb2 = -dboxAir3B[2] + fgkThPCB + fgkThPCB/2.;
  TVirtualMC::GetMC()->Gspos("EPCB", 1, "ECGB", 0., 0., zbpcb2, 0, "ONLY");
  //Positioning the Honey Comb 0.5 cm
  zhc2 = zbpcb2 + fgkThPCB/2. + fgkCellDepth/2.;
  TVirtualMC::GetMC()->Gspos("EHC2", 1, "ECGB", 0., 0., zhc2, 0, "ONLY");
  //Positioning the Top PCB 0.16 cm
  ztpcb2 = zhc2 + fgkCellDepth/2 + fgkThPCB/2.;
  TVirtualMC::GetMC()->Gspos("EPCB", 2, "ECGB", 0., 0., ztpcb2, 0, "ONLY");



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
  TVirtualMC::GetMC()->Gsvolu("EUM2","BOX", idtmed[698], dboxUM2, 3);

  //VETO
  TVirtualMC::GetMC()->Gsvolu("EUV2","BOX", idtmed[698], dboxUM2, 3);
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
  TVirtualMC::GetMC()->Gsvolu("EBPB","BOX", idtmed[607], dboxBaseB, 3);
  //----------------------------------------------------//  

  //VETO
  //- Placing of all components of UM in AIR BOX EUM2--//
  //(1)   FIRST PUT THE BASE PLATE
  Float_t zbaseplate2 = -dboxUM2[2] + fgkThBase/2.;
  TVirtualMC::GetMC()->Gspos("EBPB", 1, "EUV2", 0., 0., zbaseplate2, 0, "ONLY");

  //(2)   NEXT PLACING the SS BOX 
  Float_t zss2 = zbaseplate2 + fgkThBase/2. + fgkThSS/2.;
  TVirtualMC::GetMC()->Gspos("ESVB", 1, "EUV2", 0., 0., zss2, 0, "ONLY");
  
  // (3) Positioning the Backplane PCB 0.1 cm
  Float_t zbkp2 = zss2 + fgkThSS/2. + fgkThBKP/2.;
  TVirtualMC::GetMC()->Gspos("EBKB", 1, "EUV2", 0., 0., zbkp2, 0, "ONLY");



  //FOR PRESHOWER
  // (3) Positioning the Backplane PCB 0.1 cm
  zbkp2 = -dboxUM2[2] + fgkThBKP/2.;
  TVirtualMC::GetMC()->Gspos("EBKB", 1, "EUM2", 0., 0., zbkp2, 0, "ONLY");

  //(2)   NEXT PLACING the SS BOX 
  zss2 = zbkp2 + fgkThBKP/2. + fgkThSS/2.;
  TVirtualMC::GetMC()->Gspos("ESSB", 1, "EUM2", 0., 0., zss2, 0, "ONLY");
  
  //(1)   FIRST PUT THE BASE PLATE
  zbaseplate2 = zss2 + fgkThSS/2 + fgkThBase/2.;
  TVirtualMC::GetMC()->Gspos("EBPB", 1, "EUM2", 0., 0., zbaseplate2, 0, "ONLY");
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
  TVirtualMC::GetMC()->Gsvolu("EPB1","BOX", idtmed[600], dboxPba, 3);

  Float_t dboxPbb[3];
  dboxPbb[0] = dboxUM2[0];
  dboxPbb[1] = dboxUM2[1];
  dboxPbb[2] = fgkThLead/2.;
  // Lead of UM dimension
  TVirtualMC::GetMC()->Gsvolu("EPB2","BOX", idtmed[600], dboxPbb, 3);

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
  TVirtualMC::GetMC()->Gsvolu("ESMA","BOX", idtmed[698], dboxSM1, 3);
  
  //FOR VETO
  TVirtualMC::GetMC()->Gsvolu("EMVA","BOX", idtmed[698], dboxSM1, 3);

  //Position the 6 unit modules in EMSA
  Float_t xa1,xa2,xa3,ya1,ya2; 
  xa1 =  dboxSM1[0] - dboxUM1[0];
  xa2 = xa1 - dboxUM1[0] - 0.15 - dboxUM1[0];
  xa3 = xa2 - dboxUM1[0] - 0.15 - dboxUM1[0];
  ya1 = dboxSM1[1]  - dboxUM1[1];
  ya2 = ya1 - dboxUM1[1] - 0.1 - dboxUM1[1];

  //PRESHOWER
  // TVirtualMC::GetMC()->Gspos("EUM1", 1, "ESMA", xa1, ya1, 0., 0, "ONLY"); // BKN
  TVirtualMC::GetMC()->Gspos("EUM1", 2, "ESMA", xa2, ya1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUM1", 3, "ESMA", xa3, ya1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUM1", 4, "ESMA", xa1, ya2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUM1", 5, "ESMA", xa2, ya2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUM1", 6, "ESMA", xa3, ya2, 0., 0, "ONLY");

  //VETO
  TVirtualMC::GetMC()->Gspos("EUV1", 1, "EMVA", xa1, ya1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUV1", 2, "EMVA", xa2, ya1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUV1", 3, "EMVA", xa3, ya1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUV1", 4, "EMVA", xa1, ya2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUV1", 5, "EMVA", xa2, ya2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUV1", 6, "EMVA", xa3, ya2, 0., 0, "ONLY");


  // volume for SUPERMODULE ESMB 
  //Space is added to provide a gapping for HV between UM's
  Float_t dboxSM2[3];
  dboxSM2[0] = 2.0*dboxUM2[0] + 0.075; 
  dboxSM2[1] = 3.0*dboxUM2[1] + (2.0*0.05);
  dboxSM2[2] = dboxUM2[2];
  
  //PRESHOWER
  TVirtualMC::GetMC()->Gsvolu("ESMB","BOX", idtmed[698], dboxSM2, 3);
  //VETO 
  TVirtualMC::GetMC()->Gsvolu("EMVB","BOX", idtmed[698], dboxSM2, 3);

  //Position the 6 unit modules in EMSB
  Float_t xb1,xb2,yb1,yb2,yb3; 
  xb1 = dboxSM2[0] - dboxUM2[0];
  xb2 = xb1 - dboxUM2[0] - 0.15 - dboxUM2[0];
  yb1 = dboxSM2[1] - dboxUM2[1];
  yb2 = yb1 - dboxUM2[1] - 0.1 -  dboxUM2[1];
  yb3 = yb2 - dboxUM2[1] - 0.1 -  dboxUM2[1];


  //PRESHOWER  
  // TVirtualMC::GetMC()->Gspos("EUM2", 1, "ESMB", xb1, yb1, 0., 0, "ONLY");  // BKN
  // TVirtualMC::GetMC()->Gspos("EUM2", 2, "ESMB", xb2, yb1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUM2", 3, "ESMB", xb1, yb2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUM2", 4, "ESMB", xb2, yb2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUM2", 5, "ESMB", xb1, yb3, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUM2", 6, "ESMB", xb2, yb3, 0., 0, "ONLY");
  
  //VETO
  TVirtualMC::GetMC()->Gspos("EUV2", 1, "EMVB", xb1, yb1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUV2", 2, "EMVB", xb2, yb1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUV2", 3, "EMVB", xb1, yb2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUV2", 4, "EMVB", xb2, yb2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUV2", 5, "EMVB", xb1, yb3, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EUV2", 6, "EMVB", xb2, yb3, 0., 0, "ONLY");
  
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
  
  TVirtualMC::GetMC()->Gsvolu("ESPA","BOX", idtmed[698], dboxSMPb1, 3);
  

  //Position the 6 unit modules in ESMPbA
  Float_t xpa1,xpa2,xpa3,ypa1,ypa2; 
  xpa1 = -dboxSMPb1[0] + dboxUM1[0];
  xpa2 = xpa1 + dboxUM1[0] + 0.15 + dboxUM1[0];
  xpa3 = xpa2 + dboxUM1[0] + 0.15 + dboxUM1[0];
  ypa1 = dboxSMPb1[1]  - dboxUM1[1];
  ypa2 = ypa1 - dboxUM1[1] - 0.1 - dboxUM1[1];


  TVirtualMC::GetMC()->Gspos("EPB1", 1, "ESPA", xpa1, ypa1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPB1", 2, "ESPA", xpa2, ypa1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPB1", 3, "ESPA", xpa3, ypa1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPB1", 4, "ESPA", xpa1, ypa2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPB1", 5, "ESPA", xpa2, ypa2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPB1", 6, "ESPA", xpa3, ypa2, 0., 0, "ONLY");


  // volume for SUPERMODULE ESMPbB 
  //Space is added to provide a gapping for HV between UM's
  Float_t dboxSMPb2[3];
  dboxSMPb2[0] = 2.0*dboxUM2[0] + 0.075;
  dboxSMPb2[1] = 3.0*dboxUM2[1] + (2.0*0.05);
  dboxSMPb2[2] = fgkThLead/2.;

  TVirtualMC::GetMC()->Gsvolu("ESPB","BOX", idtmed[698], dboxSMPb2, 3);
 
  //Position the 6 unit modules in ESMPbB
  Float_t xpb1,xpb2,ypb1,ypb2,ypb3; 
  xpb1 = -dboxSMPb2[0] + dboxUM2[0];
  xpb2 = xpb1 + dboxUM2[0] + 0.15 + dboxUM2[0];
  ypb1 = dboxSMPb2[1]  - dboxUM2[1];
  ypb2 = ypb1 - dboxUM2[1] - 0.1 -  dboxUM2[1];
  ypb3 = ypb2 - dboxUM2[1] - 0.1 -  dboxUM2[1];


  TVirtualMC::GetMC()->Gspos("EPB2", 1, "ESPB", xpb1, ypb1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPB2", 2, "ESPB", xpb2, ypb1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPB2", 3, "ESPB", xpb1, ypb2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPB2", 4, "ESPB", xpb2, ypb2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPB2", 5, "ESPB", xpb1, ypb3, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPB2", 6, "ESPB", xpb2, ypb3, 0., 0, "ONLY");


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

  TVirtualMC::GetMC()->Gsvolu("EFEE","BOX", idtmed[607], dboxFEE, 3);

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
  TVirtualMC::GetMC()->Gsvolu("EFBA","BOX", idtmed[698], dboxFEEBPlaneA, 3);

  //TYPE B
  Float_t dboxFEEBPlaneB[3];
  dboxFEEBPlaneB[0]   = dboxBPlaneB[0]; //dbox4[0]; 
  dboxFEEBPlaneB[1]   = dboxBPlaneB[1];//dbox4[1];       
  dboxFEEBPlaneB[2]   = 1.2;
  
  //Back PLane PCB of MAterial G10
  TVirtualMC::GetMC()->Gsvolu("EFBB","BOX", idtmed[698], dboxFEEBPlaneB, 3);

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
	  TVirtualMC::GetMC()->Gspos("EFEE", number, "EFBA", xFee,yFee,zFee, 0, "ONLY");
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
	  TVirtualMC::GetMC()->Gspos("EFEE", number, "EFBB", xFee,yFee,zFee, 0, "ONLY");
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
  TVirtualMC::GetMC()->Gsvolu("EFSA","BOX", idtmed[698],dboxEFSA, 3);

  //Distance between the two backplanes of two UMs
  //in x-direction is 0.92 and ydirection is 0.95
  Float_t dboxEFSB[3];
  dboxEFSB[0] = 2.0*dboxFEEBPlaneB[0] + (0.938/2.0);
  dboxEFSB[1] = 3.0*dboxFEEBPlaneB[1] + 1.05;
  dboxEFSB[2] = dboxFEEBPlaneB[2];

  //Type A
  TVirtualMC::GetMC()->Gsvolu("EFSB","BOX", idtmed[698],dboxEFSB, 3);


  Float_t xfs1,xfs2,xfs3,yfs1,yfs2,yfs3; 
  xfs1 = -dboxEFSA[0] + dboxFEEBPlaneA[0];
  xfs2 = xfs1 + dboxFEEBPlaneA[0] +  0.92 + dboxFEEBPlaneA[0];
  xfs3 = xfs2 + dboxFEEBPlaneA[0] +  0.92 + dboxFEEBPlaneA[0];
  yfs1 = dboxEFSA[1] - dboxFEEBPlaneA[1];
  yfs2 = yfs1 - dboxFEEBPlaneA[1] - 0.95 - dboxFEEBPlaneA[1];



  // TVirtualMC::GetMC()->Gspos("EFBA", 1, "EFSA", xfs1, yfs1, 0., 0, "ONLY");  // BKN
  TVirtualMC::GetMC()->Gspos("EFBA", 2, "EFSA", xfs2, yfs1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EFBA", 3, "EFSA", xfs3, yfs1, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EFBA", 4, "EFSA", xfs1, yfs2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EFBA", 5, "EFSA", xfs2, yfs2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EFBA", 6, "EFSA", xfs3, yfs2, 0., 0, "ONLY");


  //Type B positioning

  xfs1 = -dboxEFSB[0] + dboxFEEBPlaneB[0];
  xfs2 = xfs1 + dboxFEEBPlaneB[0] + 0.938 + dboxFEEBPlaneB[0];
  yfs1 = dboxEFSB[1] - dboxFEEBPlaneB[1];
  yfs2 = yfs1 - dboxFEEBPlaneB[1] - 1.05 - dboxFEEBPlaneB[1];
  yfs3 = yfs2 - dboxFEEBPlaneB[1] - 1.05 - dboxFEEBPlaneB[1];



  // TVirtualMC::GetMC()->Gspos("EFBB", 1, "EFSB", xfs1, yfs1, 0., 0, "ONLY"); // BKN
  // TVirtualMC::GetMC()->Gspos("EFBB", 2, "EFSB", xfs2, yfs1, 0., 0, "ONLY"); // BKN
  TVirtualMC::GetMC()->Gspos("EFBB", 3, "EFSB", xfs1, yfs2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EFBB", 4, "EFSB", xfs2, yfs2, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EFBB", 5, "EFSB", xfs1, yfs3, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EFBB", 6, "EFSB", xfs2, yfs3, 0., 0, "ONLY");


}
 
//_____________________________________________________________________________

void AliPMDv2008::CreatePMD()
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
  
  TVirtualMC::GetMC()->Gsvolu("EFEA","BOX", idtmed[618], dboxFea, 3);

  // --- DEFINE Iron volumes  for SM B
  
  //   Fe Support 
  Float_t dboxFeb[3];
  dboxFeb[0] = fSMLengthbx;
  dboxFeb[1] = fSMLengthby;
  dboxFeb[2] = fgkThSteel/2.;
  
  TVirtualMC::GetMC()->Gsvolu("EFEB","BOX", idtmed[618], dboxFeb, 3);

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

  TVirtualMC::GetMC()->Gsvolu("EPM1", "BOX", idtmed[698], gaspmd, 3);
  TVirtualMC::GetMC()->Gsvolu("EPM2", "BOX", idtmed[698], gaspmd, 3);

  //Complete detector for Type A
  //Position Super modules type A for both CPV and PMD in EPMD  
  Float_t zpsa,zpba,zfea,zcva,zfee; 

  // zpsa = - gaspmd[2] + fSMthick/2.;
  // -2.5 is given to place PMD at -361.5 
  // BM : In future after putting proper electronics
  // -2.5 will be replaced by -gaspmd[2]

  //TYPE A
  //Fee board

  // This part is commented for the time being by BKN

  zfee=-gaspmd[2] + 1.2;

  /*
  TVirtualMC::GetMC()->Gspos("EFSA", 1, "EPM1", 0., 0., zfee, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EFSA", 2, "EPM2", 0., 0., zfee, jhrot12, "ONLY");
  */

  //VETO

  zcva = zfee + 1.2 + fDthick;

  /*
  TVirtualMC::GetMC()->Gspos("EMVA", 1, "EPM1", 0., 0., zcva, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EMVA", 2, "EPM2", 0., 0., zcva, jhrot12, "ONLY");
  */



  //Iron support
  zfea = zcva + fDthick + fgkThSteel/2.;
  TVirtualMC::GetMC()->Gspos("EFEA", 1, "EPM1", 0., 0., zfea, 0, "ONLY");
  //TVirtualMC::GetMC()->Gspos("EFEA", 2, "EPM2", 0., 0., zfea, 0, "ONLY");
  //Lead
  zpba=zfea+fgkThSteel/2.+ fgkThLead/2.;
  TVirtualMC::GetMC()->Gspos("ESPA", 1, "EPM1", 0., 0., zpba, 0, "ONLY");
  //TVirtualMC::GetMC()->Gspos("ESPA", 2, "EPM2", 0., 0., zpba, 0, "ONLY");
  //Preshower
  zpsa = zpba + fgkThLead/2. + fDthick;
  TVirtualMC::GetMC()->Gspos("ESMA", 1, "EPM1", 0., 0., zpsa, 0, "ONLY");
  //TVirtualMC::GetMC()->Gspos("ESMA", 2, "EPM2", 0., 0., zpsa, jhrot12, "ONLY");
  //FEE boards
  zfee=zpsa + fDthick + 1.2;
  TVirtualMC::GetMC()->Gspos("EFSA", 3, "EPM1", 0., 0., zfee, 0, "ONLY");
  //TVirtualMC::GetMC()->Gspos("EFSA", 4, "EPM2", 0., 0., zfee, jhrot12, "ONLY");

 
  //TYPE - B
  gaspmd[0] = fSMLengthbx; 
  gaspmd[1] = fSMLengthby; 
  gaspmd[2] = fSMthick; 

  TVirtualMC::GetMC()->Gsvolu("EPM3", "BOX", idtmed[698], gaspmd, 3);
  TVirtualMC::GetMC()->Gsvolu("EPM4", "BOX", idtmed[698], gaspmd, 3);

  //Complete detector for Type B
  //Position Super modules type B for both CPV and PMD in EPMD  
  Float_t zpsb,zpbb,zfeb,zcvb; 
  // zpsb = - gaspmd[2] + fSMthick/2.;
  // -2.5 is given to place PMD at -361.5 
  // BM: In future after putting proper electronics
  // -2.5 will be replaced by -gaspmd[2]

 //Fee board

  zfee=-gaspmd[2] + 1.2;

  /*
  TVirtualMC::GetMC()->Gspos("EFSB", 5, "EPM3", 0., 0., zfee, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EFSB", 6, "EPM4", 0., 0., zfee, jhrot12, "ONLY");
  */

  zcvb= zfee + 1.2 + fDthick;

  //VETO
  /*
  TVirtualMC::GetMC()->Gspos("EMVB", 3, "EPM3", 0., 0., zcvb, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EMVB", 4, "EPM4", 0., 0., zcvb, jhrot12, "ONLY");
  */

  //IRON SUPPORT
  zfeb= zcvb + fDthick +  fgkThSteel/2.;
  //TVirtualMC::GetMC()->Gspos("EFEB", 3, "EPM3", 0., 0., zfeb, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EFEB", 4, "EPM4", 0., 0., zfeb, 0, "ONLY");
  //LEAD
  zpbb= zfeb + fgkThSteel/2.+ fgkThLead/2.;
  //TVirtualMC::GetMC()->Gspos("ESPB", 3, "EPM3", 0., 0., zpbb, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ESPB", 4, "EPM4", 0., 0., zpbb, 0, "ONLY");
  //PRESHOWER
  zpsb = zpbb + fgkThLead/2.+ fDthick;
  //TVirtualMC::GetMC()->Gspos("ESMB", 3, "EPM3", 0., 0., zpsb, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ESMB", 4, "EPM4", 0., 0., zpsb, jhrot12, "ONLY");
  //FEE boards
  zfee=zpsb + fDthick + 1.2;
  //TVirtualMC::GetMC()->Gspos("EFSB", 7, "EPM3", 0., 0., zfee, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EFSB", 8, "EPM4", 0., 0., zfee, jhrot12, "ONLY");


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
  TVirtualMC::GetMC()->Gspos("EPM1", 1, "ALIC",  xsma,ysma,zp,  0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPM2", 1, "ALIC", -xsma,-ysma,zp, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPM3", 1, "ALIC",  xsmb,ysmb,zp,  0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPM4", 1, "ALIC", -xsmb,-ysmb,zp, 0, "ONLY");
}

 
//_____________________________________________________________________________
void AliPMDv2008::CreateMaterials()
{
  // Create materials for the PMD
  //
  // ORIGIN    : Y. P. VIYOGI 
  //
  //  cout << " Inside create materials " << endl;

  Int_t isxfld = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
  Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();
  
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
  
  AliDebug(1,"Outside create materials");

}

//_____________________________________________________________________________
void AliPMDv2008::Init()
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
  // --- Generate explicitly delta rays in the iron, aluminium and lead --- 
  // Gstpar removed from here and all energy cut-offs moved to galice.cuts
  // Visualization of volumes
  gGeoManager->SetVolumeAttribute("ECAR", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECCU", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECCU", "COLO", 4);
  gGeoManager->SetVolumeAttribute("EST1", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EST2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EHC1", "SEEN", 0);  
  gGeoManager->SetVolumeAttribute("EHC2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EPCA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EBKA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECGA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECVA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EDGA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EDVA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ESSA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ESVA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EUM1", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EUV1", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EBPA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EPCB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EBKB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECGB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECVB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EDGB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EDVB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ESSB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ESVB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EUM2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EUV2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EBPB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EPB1", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EPB2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ESMA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EMVA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ESMB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EMVB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ESPA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ESPB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EFEE", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EFEE", "COLO", 4);
  gGeoManager->SetVolumeAttribute("EFBA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EFBB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EFSA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EFSB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EFEA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EFEB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EPM1", "SEEN", 1);
  gGeoManager->SetVolumeAttribute("EPM2", "SEEN", 1);
  gGeoManager->SetVolumeAttribute("EPM3", "SEEN", 1);
  gGeoManager->SetVolumeAttribute("EPM4", "SEEN", 1);
}

//_____________________________________________________________________________
void AliPMDv2008::StepManager()
{
  //
  // Called at each step in the PMD
  //

  Int_t   copy;
  Float_t hits[5], destep;
  Float_t center[3] = {0,0,0};
  Int_t   vol[6];
  
  if(fMC->CurrentMedium() == fMedSens && (destep = fMC->Edep())) {
  
    fMC->CurrentVolID(copy);
    vol[0] = copy;

    fMC->CurrentVolOffID(1,copy);
    vol[1] = copy;

    fMC->CurrentVolOffID(2,copy);
    vol[2] = copy;

    fMC->CurrentVolOffID(3,copy);
    vol[3] = copy;

    fMC->CurrentVolOffID(4,copy);
    vol[4] = copy;

    fMC->CurrentVolOffID(5,copy);
    vol[5] = copy;


    fMC->Gdtom(center,hits,1);
    hits[3] = destep*1e9; //Number in eV

    // this is for pile-up events
    hits[4] = fMC->TrackTime();

    AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);

  }
}

  
//------------------------------------------------------------------------
// Get parameters

void AliPMDv2008::GetParameters()
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
void AliPMDv2008::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  // 
  SetSectorAlignable();

}
// ----------------------------------------------------------------
void AliPMDv2008::SetSectorAlignable() const
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
    if(!gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data()))
      {
	AliFatal("Unable to set alignable entry!");
      }
  }
}
// ------------------------------------------------------------------
