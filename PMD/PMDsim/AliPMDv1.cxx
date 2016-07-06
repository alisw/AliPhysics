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
//---------------------------------------------------     
//  ALICE PMD FEE BOARDS IMPLEMENTATION
//  Dt: 25th February 2006 
//  M.M. Mondal, S.K. Prasad and P.K. Netrakanti
//---------------------------------------------------
//   Create final detector from Unit Modules
//   Author : Bedanga and Viyogi June 2003
//---------------------------------------------------
// Modified by
// Dr. Y.P. Viyogi and Ranbir Singh
// Dt: 2nd February 2009
//
//Begin_Html
/*
<img src="picts/AliPMDv1Class.gif">
*/
//End_Html
//                                                                           //
/////////////////////////////////////////////////////////////////////////////
////

#include <Riostream.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>
#include <TVirtualMC.h>

#include "AliConst.h" 
#include "AliLog.h"
#include "AliMC.h"
#include "AliMagF.h" 
#include "AliPMDv1.h"
#include "AliRun.h"
#include "AliTrackReference.h"

const Int_t   AliPMDv1::fgkNcolUM1    = 48;     // Number of cols in UM, type 1
const Int_t   AliPMDv1::fgkNcolUM2    = 96;     // Number of cols in UM, type 2
const Int_t   AliPMDv1::fgkNrowUM1    = 96;     // Number of rows in UM, type 1
const Int_t   AliPMDv1::fgkNrowUM2    = 48;     // Number of rows in UM, type 2
const Float_t AliPMDv1::fgkCellRadius = 0.25;     // Radius of a hexagonal cell
const Float_t AliPMDv1::fgkCellWall   = 0.02;     // Thickness of cell Wall
const Float_t AliPMDv1::fgkCellDepth  = 0.50;     // Gas thickness
const Float_t AliPMDv1::fgkThPCB      = 0.16;     // Thickness of PCB 
const Float_t AliPMDv1::fgkThLead     = 1.5;      // Thickness of Pb
const Float_t AliPMDv1::fgkThSteel    = 0.5;      // Thickness of Steel
const Float_t AliPMDv1::fgkGap        = 0.025;    // Air Gap
const Float_t AliPMDv1::fgkZdist      = 361.5;    // z-position of the detector
const Float_t AliPMDv1::fgkSqroot3    = 1.7320508;// Square Root of 3
const Float_t AliPMDv1::fgkSqroot3by2 = 0.8660254;// Square Root of 3 by 2
const Float_t AliPMDv1::fgkSSBoundary = 0.3;
const Float_t AliPMDv1::fgkThSS       = 1.23;     // Old thickness of SS frame was 1.03
const Float_t AliPMDv1::fgkThTopG10   = 0.33;
const Float_t AliPMDv1::fgkThBotG10   = 0.4;


ClassImp(AliPMDv1)
 
//_____________________________________________________________________________
AliPMDv1::AliPMDv1():
  fSMthick(0.),
  fSMthickpmd(0.),
  fDthick(0.),
  fSMLengthax(0.),
  fSMLengthay(0.),
  fSMLengthbx(0.),
  fSMLengthby(0.),
  fMedSens(0)
{
  
  // Default constructor 
  
  for (Int_t i = 0; i < 3; i++)
    {
      fDboxmm1[i]  = 0.;
      fDboxmm12[i] = 0.;
      fDboxmm2[i]  = 0.;
      fDboxmm22[i] = 0.;
    }
  for (Int_t i = 0; i < 48; i++)
    {
      fModStatus[i] = 1;
    }

}
 
//_____________________________________________________________________________
AliPMDv1::AliPMDv1(const char *name, const char *title):
  AliPMD(name,title),
  fSMthick(0.),
  fSMthickpmd(0.),
  fDthick(0.),
  fSMLengthax(0.),
  fSMLengthay(0.),
  fSMLengthbx(0.),
  fSMLengthby(0.),
  fMedSens(0)
{
  
  // Standard constructor
  
  for (Int_t i = 0; i < 3; i++)
    {
      fDboxmm1[i]  = 0.;
      fDboxmm12[i] = 0.;
      fDboxmm2[i]  = 0.;
      fDboxmm22[i] = 0.;
    }
  for (Int_t i = 0; i < 48; i++)
    {
      fModStatus[i] = 1;
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
  // Creates the geometry of the cells of PMD, places them in  modules 
  // which are rectangular objects.
  // Basic unit is ECAR, a hexagonal cell made of Ar+CO2, which is 
  // placed inside another hexagonal cell made of Cu (ECCU) with larger 
  // radius, compared to ECAR. The difference in radius gives the dimension 
  // of half width of each cell wall.
  // These cells are placed in a rectangular strip which are of 2 types 
  // EST1 and EST2. 
  // Two types of honeycomb EHC1 & EHC2 are made using strips EST1 & EST2. 
  // 4 types of unit modules are made EUM1 & EUM2 for PRESHOWER Plane and
  // EUV1 & EUV2 for VETO Plane which contains  strips placed repeatedly 
  //  
  // These unit moules are then placed inside EPM1, EPM2, EPM3 and EPM4 along
  // with lead convertor ELDA & ELDB and Iron Supports EFE1, EFE2, EFE3 and EFE4
  //  They have 6 unit moudles inside them in each plane. Therefore, total of 48
  // unit modules in both the planes (PRESHOWER Plane & VETO Plane). The numbering
  // of unit modules is from 0 to 47.
  //
  // Steel channels (ECHA & ECHB) are also placed which are used to place the unit modules
  // 
  // In order to account for the extra material around and on the detector, Girders (EGDR),
  // girder's Carriage (EXGD), eight Aluminium boxes (ESV1,2,3,4 & EVV1,2,3,4) along with
  // LVDBs (ELVD), cables (ECB1,2,3,4), and ELMBs (ELMB) are being placed in approximations.
  // 
  //  Four FR4 sheets (ECC1,2,3,4) are placed parallel to the PMD on both sides, which perform 
  // as cooling encloser
 
  // NOTE:-  VOLUME Names : begining with "E" for all PMD volumes 
  
  Int_t i,j;
  Int_t number;
  Int_t ihrotm,irotdm;
  Float_t xb, yb, zb;

  Int_t *idtmed = fIdtmed->GetArray()-599;
 
  AliMatrix(ihrotm, 90., 30.,   90.,  120., 0., 0.);
  AliMatrix(irotdm, 90., 180.,  90.,  270., 180., 0.);
 
  //******************************************************//
  //                    STEP - I                          //
  //******************************************************//
  // First create the sensitive medium of a hexagon cell (ECAR)
  // Inner hexagon filled with gas (Ar+CO2)
  // Integer assigned to Ar+CO2 medium is 604

  Float_t hexd2[10] = {0.,360.,6,2,-0.25,0.,0.23,0.25,0.,0.23};
  hexd2[4] = -fgkCellDepth/2.;
  hexd2[7] =  fgkCellDepth/2.;
  hexd2[6] =  fgkCellRadius - fgkCellWall;
  hexd2[9] =  fgkCellRadius - fgkCellWall;
  
  TVirtualMC::GetMC()->Gsvolu("ECAR", "PGON", idtmed[604], hexd2,10);

  //******************************************************//
  //                    STEP - II                         //
  //******************************************************//
  // Place the sensitive medium inside a hexagon copper cell (ECCU)
  // Outer hexagon made of Copper
  // Integer assigned to Cu medium is 614
  
  Float_t hexd1[10] = {0.,360.,6,2,-0.25,0.,0.25,0.25,0.,0.25};
  hexd1[4] = -fgkCellDepth/2.;
  hexd1[7] =  fgkCellDepth/2.;
  hexd1[6] =  fgkCellRadius;
  hexd1[9] =  fgkCellRadius;
  
  TVirtualMC::GetMC()->Gsvolu("ECCU", "PGON", idtmed[614], hexd1,10);

  // Place  inner hex (sensitive volume) inside outer hex (copper)
  
  TVirtualMC::GetMC()->Gspos("ECAR", 1, "ECCU", 0., 0., 0., 0, "ONLY");

  //******************************************************//
  //                    STEP - III                        //
  //******************************************************//
  // Now create Two types of Rectangular strips (EST1, EST2) 
  // of 1 column and 96 or 48 cells length

  // volume for first strip EST1 made of AIR 
  // Integer assigned to Air medium is 698
  // strip type-1 is of 1 column and 96 rows i.e. of 96 cells length 

  Float_t dbox1[3];
  dbox1[0] = fgkCellRadius/fgkSqroot3by2;
  dbox1[1] = fgkNrowUM1*fgkCellRadius;
  dbox1[2] = fgkCellDepth/2.;
  
  TVirtualMC::GetMC()->Gsvolu("EST1","BOX", idtmed[698], dbox1, 3);


  // volume for second strip EST2 
  // strip type-2 is of 1 column and 48 rows i.e. of 48 cells length 

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
      yb -= (fgkCellRadius*2.);
    }
  
  
  //******************************************************//
  //                    STEP - IV                         //
  //******************************************************//
  // Create EHC1 : The honey combs for a unit module type-1
  //-------------------------EHC1 Start-------------------//
  
  // First step is to create a honey comb unit module.
  // This is named as EHC1 and  is a volume of Air 
  // we will lay the EST1 strips of honey comb cells inside it.
  
  // Dimensions of EHC1
  // X-dimension = (dbox1[0]*fgkNcolUM1)-(fgkCellRadius*fgkSqroot3*(fgkNcolUM1-1)/6.)+ 0.15+0.05+0.05; 
  // Y-dimension = Number of rows * cell radius/sqrt3by2 + 0.15+0.05+0.05;  
  // 0.15cm is the extension in honeycomb on both side of X and Y, 0.05 for air gap and 0.05
  // for G10 boundary around, which are now merged in the dimensions of EHC1 
  // Z-dimension = cell depth/2

  Float_t ehcExt = 0.15;
  Float_t ehcAround = 0.05 + 0.05;;

  Float_t dbox3[3];
  dbox3[0] = (dbox1[0]*fgkNcolUM1)-
    (fgkCellRadius*fgkSqroot3*(fgkNcolUM1-1)/6.) + ehcExt + ehcAround;  
  dbox3[1] = dbox1[1]+fgkCellRadius/2. + ehcExt + ehcAround; 
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
      TVirtualMC::GetMC()->Gspos("EST1",number, "EHC1", xb - 0.25, yb , 0. , 0, "MANY");
      
      //The strips are being placed from top towards bottom of the module
      //This is because the first cell in a module in hardware is the top
      //left corner cell
      
      xb = (dbox3[0]-dbox1[0])-j*fgkCellRadius*fgkSqroot3;
      
    }
  
  //--------------------EHC1 done----------------------------------------//
  
  
  
  //--------------------------------EHC2 Start---------------------------//
  // Create EHC2 : The honey combs for a unit module type-2 
  // First step is to create a honey comb unit module.
  // This is named as EHC2, we will lay the EST2 strips of
  // honey comb cells inside it.
  
  // Dimensions of EHC2
  // X-dimension = (dbox2[0]*fgkNcolUM2)-(fgkCellRadius*fgkSqroot3*(fgkNcolUM2-1)/6.)+ 0.15+0.05+0.05;
  // Y-dimension = Number of rows * cell radius/sqrt3by2 + 0.15+0.05+0.05;
  // 0.15cm is the extension in honeycomb on both side of X and Y, 0.05 for air gap and 0.05
  // for G10 boundary around, which are now merged in the dimensions of EHC2 
  // Z-dimension = cell depth/2
  
  
  Float_t dbox4[3];
  
  dbox4[0] =(dbox2[0]*fgkNcolUM2)-
    (fgkCellRadius*fgkSqroot3*(fgkNcolUM2-1)/6.) + ehcExt + ehcAround; 
  dbox4[1] = dbox2[1] + fgkCellRadius/2. + ehcExt + ehcAround;
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
      TVirtualMC::GetMC()->Gspos("EST2",number, "EHC2", xb - 0.25, yb , 0. ,0, "MANY");
      xb = (dbox4[0]-dbox2[0])-j*fgkCellRadius*fgkSqroot3;
    }
  
 
  //----------------------------EHC2 done-------------------------------//

  //====================================================================//
 
  // Now the job is to assmeble an Unit module
  // It will have the following components
  // (a) Base plate of G10 of 0.2cm 
  // (b) Air gap  of 0.08cm   
  // (c) Bottom PCB of 0.16cm G10
  // (d) Honey comb 0f 0.5cm
  // (e) Top PCB  of 0.16cm G10 
  // (f) Back Plane of 0.1cm G10
  // (g) Then all around then we have an air gap of 0.05cm
  // (h) Then all around 0.05cm thick G10 insulation
  // (i) Then all around Stainless Steel boundary channel 0.3 cm thick

  // In order to reduce the number of volumes and simplify the geometry
  // following steps are performed:
  // (I)   Base Plate(0.2cm), Air gap(0.04cm) and Bottom PCB(0.16cm) 
  //       are taken together as a G10 Plate EDGA (0.4cm)
  // (II) Back Plane(0.1cm), Air Gap(0.04cm) and Top PCB(0.16cm) and extra 
  //      clearance 0.03cm are taken together as G10 Plate EEGA(0.33cm)
  // (III) The all around Air gap(0.05cm) and G10 boundary(0.05cm) are already 
  //       merged in the dimension of EHC1, EHC2, EDGA and EEGA. Therefore, no 
  //       separate volumes for all around materials
  
  //Let us first create them one by one
  //--------------------------------------------------------------------//

  // ---------------- Lets do it first for UM Long Type -----//
  // 4mm G10 Box : Bottom PCB + Air Gap + Base Plate
  //================================================
  // Make a 4mm thick G10 Box for Unit module Long Type 
  // X-dimension is EHC1 - ehcExt
  // Y-dimension is EHC1 - ehcExt
  // EHC1 was extended 0.15cm(ehcExt) on both sides
  // Z-dimension 0.4/2 = 0.2 cm
  // Integer assigned to G10 medium is 607

  Float_t dboxCGA[3];
  dboxCGA[0]  = dbox3[0] - ehcExt; 
  dboxCGA[1]  = dbox3[1] - ehcExt; 
  dboxCGA[2]  = fgkThBotG10/2.;

  //Create a G10 BOX 
  TVirtualMC::GetMC()->Gsvolu("EDGA","BOX", idtmed[607], dboxCGA, 3);

  //-------------------------------------------------//
  // 3.3mm G10 Box : Top PCB + Air GAp + Back Plane
  //================================================
  // Make a 3.3mm thick G10 Box for Unit module Long Type 
  // X-dimension is EHC1 - ehcExt
  // Y-dimension is EHC1 - ehcExt
  // EHC1 was extended 0.15cm(ehcExt) on both sides
  // Z-dimension 0.33/2 = 0.165 cm

  Float_t dboxEEGA[3];
  dboxEEGA[0]  = dboxCGA[0]; 
  dboxEEGA[1]  = dboxCGA[1]; 
  dboxEEGA[2]  = fgkThTopG10/2.;

  //Create a G10 BOX 
  TVirtualMC::GetMC()->Gsvolu("EEGA","BOX", idtmed[607], dboxEEGA, 3);


  //----------------------------------------------------------//
  //Stainless Steel Bounadry : EUM1 & EUV1
  //
  // Make a 3.63cm thick Stainless Steel boundary for Unit module Long Type
  // 3.63cm equivalent to EDGA(0.4cm)+EHC1(0.5cm)+EEGA(0.33cm)+FEE Board(2.4cm)
  // X-dimension is EEGA + fgkSSBoundary
  // Y-dimension is EEGA + fgkSSBoundary
  // Z-dimension 1.23/2 + 2.4/2.
  // FEE Boards are 2.4cm thick
  // Integer assigned to Stainless Steel medium is 618
  //------------------------------------------------------//
  // A Stainless Steel Boundary Channel to house the unit module
  // along with the FEE Boards

  Float_t dboxSS1[3];
  dboxSS1[0]           = dboxCGA[0]+fgkSSBoundary; 
  dboxSS1[1]           = dboxCGA[1]+fgkSSBoundary;       
  dboxSS1[2]           = fgkThSS/2.+ 2.4/2.;
  
  //FOR PRESHOWER
  //Stainless Steel boundary - Material Stainless Steel
  TVirtualMC::GetMC()->Gsvolu("EUM1","BOX", idtmed[618], dboxSS1, 3);
  
  //FOR VETO
  //Stainless Steel boundary - Material Stainless Steel
  TVirtualMC::GetMC()->Gsvolu("EUV1","BOX", idtmed[618], dboxSS1, 3);
  
  //--------------------------------------------------------------------//


  

  // ============ PMD FEE BOARDS IMPLEMENTATION ======================//
  
  // FEE board
  // It is FR4 board of length * breadth :: 7cm * 2.4 cm
  // and thickness 0.2cm
  // Material medium is same as G10

  Float_t dboxFEE[3];
  dboxFEE[0] = 0.2/2.;  
  dboxFEE[1] = 7.0/2.;
  dboxFEE[2] = 2.4/2.;

  TVirtualMC::GetMC()->Gsvolu("EFEE","BOX", idtmed[607], dboxFEE, 3);

  // Now to create the Mother volume to accomodate FEE boards
  // It should have the dimension few mm smaller than the back plane
  // But, we have taken it as big as EUM1 or EUV1
  // It is to compensate the Stainless Steel medium of EUM1 or EUV1

  // Create Mother volume of Air : Long TYPE

  Float_t dboxFEEBPlaneA[3];
  dboxFEEBPlaneA[0]   = dboxSS1[0];  
  dboxFEEBPlaneA[1]   = dboxSS1[1];
  dboxFEEBPlaneA[2]   = 2.4/2.;
  
  //Volume of same dimension as EUM1 or EUV1 of Material AIR
  TVirtualMC::GetMC()->Gsvolu("EFBA","BOX", idtmed[698], dboxFEEBPlaneA, 3);
  
  //Placing the FEE boards in the Mother volume of AIR
  

  Float_t xFee;          // X-position of FEE board
  Float_t yFee;          // Y-position of FEE board
  Float_t zFee = 0.0;    // Z-position of FEE board
  
  Float_t xA    = 0.5;   //distance from the border to 1st FEE board/Translator
  Float_t yA    = 4.00;  //distance from the border to 1st FEE board
  Float_t xSepa = 1.70;  //Distance between two FEE boards in X-side
  Float_t ySepa = 8.00;  //Distance between two FEE boards in Y-side
  
  
  
  // FEE Boards EFEE placed inside EFBA
  
  yFee =  dboxFEEBPlaneA[1] - yA - 0.1 - 0.3;
  // 0.1cm and 0.3cm are subtracted to shift the FEE Boards on their actual positions
  // As the positions are changed, because we have taken the dimension of EFBA equal 
  // to the dimension of EUM1 or EUV1  
  number = 1;
  // The loop for six rows of FEE Board
  for (i = 1; i <= 6; ++i)
    {
      // First we place the translator board
      xFee = -dboxFEEBPlaneA[0] + xA + 0.1 +0.3;
      
      TVirtualMC::GetMC()->Gspos("EFEE", number, "EFBA", xFee,yFee,zFee, 0, "ONLY");
      
      // The first FEE board is 11mm from the translator board
      xFee   += 1.1;
      number += 1;
      
      for (j = 1; j <= 12; ++j)
        {
          TVirtualMC::GetMC()->Gspos("EFEE", number, "EFBA", xFee,yFee,zFee, 0, "ONLY");
          xFee += xSepa;
          number += 1;
        }
      yFee -= ySepa;
    }
  
  
  // Now Place EEGA, EDGA, EHC1 and EFBA in EUM1 & EUV1 to complete the unit module
  
  
  //                   FOR PRE SHOWER                //
  // Placing of all components of UM in AIR BOX EUM1 //
  
  //(1)   FIRST PUT the 4mm G10 Box : EDGA
  Float_t zedga = -dboxSS1[2] + fgkThBotG10/2.;
  TVirtualMC::GetMC()->Gspos("EDGA", 1, "EUM1", 0., 0., zedga, 0, "ONLY");
  
  //(2)   NEXT PLACING the Honeycomb EHC1
  Float_t zehc1 = zedga + fgkThBotG10/2. + fgkCellDepth/2.;
  TVirtualMC::GetMC()->Gspos("EHC1", 1, "EUM1", 0., 0.,  zehc1, 0, "ONLY");
  
  //(3)   NEXT PLACING the 3.3mm G10 Box : EEGA
  Float_t zeega = zehc1 + fgkCellDepth/2. + fgkThTopG10/2.;
  TVirtualMC::GetMC()->Gspos("EEGA", 1, "EUM1", 0., 0., zeega, 0, "ONLY");
  
  //(4)   NEXT PLACING the FEE BOARD : EFBA
  Float_t zfeeboardA = zeega + fgkThTopG10/2. +1.2;
  TVirtualMC::GetMC()->Gspos("EFBA", 1, "EUM1", 0., 0., zfeeboardA, 0, "ONLY");
  
  //                    FOR VETO                       //
  //  Placing of all components of UM in AIR BOX EUV1  //
  
  //(1)  FIRST PUT the FEE BOARD : EFBA
  zfeeboardA = -dboxSS1[2] + 1.2;
  TVirtualMC::GetMC()->Gspos("EFBA", 1, "EUV1", 0., 0., zfeeboardA, 0, "ONLY");
  
  //(2)  FIRST PLACING the 3.3mm G10 Box : EEGA
  zeega = zfeeboardA + 1.2 + fgkThTopG10/2.;
  TVirtualMC::GetMC()->Gspos("EEGA", 1, "EUV1", 0., 0., zeega, 0, "ONLY");
  
  //(3)   NEXT PLACING the Honeycomb EHC1
  zehc1 = zeega + fgkThTopG10/2 + fgkCellDepth/2.;
  TVirtualMC::GetMC()->Gspos("EHC1", 1, "EUV1", 0., 0.,  zehc1, 0, "ONLY");
  
  //(4)   NEXT PUT THE 4mm G10 Box : EDGA
  zedga = zehc1 + fgkCellDepth/2.+ fgkThBotG10/2.;
  TVirtualMC::GetMC()->Gspos("EDGA", 1, "EUV1", 0., 0., zedga, 0, "ONLY");
  

  //===================  LONG TYPE COMPLETED  =========================//
  //------------ Lets do the same thing for UM Short Type -------------//
  // 4mm G10 Box : Bottom PCB + Air Gap + Base Plate
  //================================================
  // Make a 4mm thick G10 Box for Unit module ShortType
  // X-dimension is EHC2 - ehcExt
  // Y-dimension is EHC2 - ehcExt
  // EHC2 was extended 0.15cm(ehcExt) on both sides
  // Z-dimension 0.4/2 = 0.2 cm
  // Integer assigned to G10 medium is 607
  
  Float_t dboxCGB[3];
  dboxCGB[0]  = dbox4[0] - ehcExt; 
  dboxCGB[1]  = dbox4[1] - ehcExt; 
  dboxCGB[2]  = 0.4/2.;
  
  //Create a G10 BOX 
  TVirtualMC::GetMC()->Gsvolu("EDGB","BOX", idtmed[607], dboxCGB, 3);
  
  //-------------------------------------------------//
  // 3.3mm G10 Box : PCB + Air Gap + Back Plane
  //================================================
  // Make a 3.3mm thick G10 Box for Unit module Short Type 
  // X-dimension is EHC2 - ehcExt
  // Y-dimension is EHC2 - ehcExt
  // EHC2 was extended 0.15cm(ehcExt) on both sides
  // Z-dimension 0.33/2 = 0.165 cm
  
  Float_t dboxEEGB[3];
  dboxEEGB[0]  = dboxCGB[0]; 
  dboxEEGB[1]  = dboxCGB[1]; 
  dboxEEGB[2]  = 0.33/2.;
  
  // Create a G10 BOX 
  TVirtualMC::GetMC()->Gsvolu("EEGB","BOX", idtmed[607], dboxEEGB, 3);
  
  
  //Stainless Steel Bounadry : EUM2 & EUV2
  //==================================
  // Make a 3.63cm thick Stainless Steel boundary for Unit module Short Type 
  // 3.63cm equivalent to EDGB(0.4cm)+EHC2(0.5cm)+EEGB(0.33cm)+FEE Board(2.4cm)
  // X-dimension is EEGB + fgkSSBoundary
  // Y-dimension is EEGB + fgkSSBoundary
  // Z-dimension 1.23/2 + 2.4/2.
  // FEE Boards are 2.4cm thick
  // Integer assigned to Stainless Steel medium is 618
  //------------------------------------------------------//
  // A Stainless Steel Boundary Channel to house the unit module
  // along with the FEE Boards
  
  
  Float_t dboxSS2[3];
  dboxSS2[0]  = dboxCGB[0] + fgkSSBoundary; 
  dboxSS2[1]  = dboxCGB[1] + fgkSSBoundary;       
  dboxSS2[2]  = fgkThSS/2.+ 2.4/2.;
  
  //PRESHOWER
  //Stainless Steel boundary - Material Stainless Steel
  TVirtualMC::GetMC()->Gsvolu("EUM2","BOX", idtmed[618], dboxSS2, 3);
  
  //VETO
  //Stainless Steel boundary - Material Stainless Steel
  TVirtualMC::GetMC()->Gsvolu("EUV2","BOX", idtmed[618], dboxSS2, 3);
  
  //----------------------------------------------------------------//
  //NOW THE FEE BOARD IMPLEMENTATION
  
  // To create the Mother volume to accomodate FEE boards
  // It should have the dimension few mm smaller than the back plane
  // But, we have taken it as big as EUM2 or EUV2
  // It is to compensate the Stainless Steel medium of EUM2 or EUV2

  // Create Mother volume of Air : SHORT TYPE 
  //------------------------------------------------------//


  Float_t dboxFEEBPlaneB[3];
  dboxFEEBPlaneB[0]   = dboxSS2[0];  
  dboxFEEBPlaneB[1]   = dboxSS2[1];       
  dboxFEEBPlaneB[2]   = 2.4/2.;
  
  //Volume of same dimension as EUM2 or EUV2 of Material AIR
  TVirtualMC::GetMC()->Gsvolu("EFBB","BOX", idtmed[698], dboxFEEBPlaneB, 3);
  
  
  // FEE Boards EFEE placed inside EFBB
  
  yFee =  dboxFEEBPlaneB[1] - yA -0.1 -0.3;  
  // 0.1cm and 0.3cm are subtracted to shift the FEE Boards on their actual positions
  // As the positions are changed, because we have taken the dimension of EFBB equal 
  // to the dimension of EUM2 or EUV2  
  number = 1;
  for (i = 1; i <= 3; ++i) 
    {
      xFee = -dboxFEEBPlaneB[0] + xA + 0.1 +0.3;  
      
      //First we place the translator board
      TVirtualMC::GetMC()->Gspos("EFEE", number, "EFBB", xFee,yFee,zFee, 0, "ONLY");
      // The first FEE board is 11mm from the translator board    
      xFee+=1.1;
      number+=1;
      
      for (j = 1; j <= 12; ++j) 
	{
	  TVirtualMC::GetMC()->Gspos("EFEE", number, "EFBB", xFee,yFee,zFee, 0, "ONLY");
	  xFee += xSepa;
	  number += 1;
	}
      
      //Now we place Bridge Board
      xFee = xFee - xSepa + 0.8 ;
      //Bridge Board is at a distance 8mm from FEE board
      TVirtualMC::GetMC()->Gspos("EFEE", number, "EFBB", xFee,yFee,zFee, 0, "ONLY");
      
      number+=1;
      xFee+=0.8;
      
      for (j = 1; j <= 12; ++j) 
	{
	  TVirtualMC::GetMC()->Gspos("EFEE", number, "EFBB", xFee,yFee,zFee, 0, "ONLY");
	  xFee += xSepa;
	  number += 1;
	}
      yFee -= ySepa; 
    }
  
  
  
  // Now Place EEGB, EDGB, EHC2 and EFBB in EUM2 & EUV2 to complete the unit module
  
  // FOR PRE SHOWER
  //- Placing of all components of UM in AIR BOX EUM2--//
  //(1)   FIRST PUT the G10 Box : EDGB
  Float_t zedgb = -dboxSS2[2] + 0.4/2.;
  TVirtualMC::GetMC()->Gspos("EDGB", 1, "EUM2", 0., 0., zedgb, 0, "ONLY");
  
  //(2)   NEXT PLACING the Honeycomb EHC2
  Float_t zehc2 = zedgb + 0.4/2. + fgkCellDepth/2.;
  TVirtualMC::GetMC()->Gspos("EHC2", 1, "EUM2", 0., 0.,  zehc2, 0, "ONLY");
  
  //(3)   NEXT PLACING the G10 Box : EEGB
  Float_t zeegb = zehc2 + fgkCellDepth/2. + 0.33/2.;
  TVirtualMC::GetMC()->Gspos("EEGB", 1, "EUM2", 0., 0., zeegb, 0, "ONLY");
  
  //(4)   NEXT PLACING FEE BOARDS : EFBB
  Float_t zfeeboardB = zeegb + 0.33/2.+1.2;
  TVirtualMC::GetMC()->Gspos("EFBB", 1, "EUM2", 0., 0., zfeeboardB, 0, "ONLY");
  
  //  FOR VETO
  //  Placing of all components of UM in AIR BOX EUV2 //
  
  //(1)  FIRST PUT the FEE BOARD : EUV2
  zfeeboardB = -dboxSS2[2] + 1.2;
  TVirtualMC::GetMC()->Gspos("EFBB", 1, "EUV2", 0., 0., zfeeboardB, 0, "ONLY");
  
  //(2)  FIRST PLACING the G10 Box : EEGB
  zeegb = zfeeboardB + 1.2 + 0.33/2.;
  TVirtualMC::GetMC()->Gspos("EEGB", 1, "EUV2", 0., 0., zeegb, 0, "ONLY");
  
  //(3)   NEXT PLACING the Honeycomb EHC2
  zehc2 = zeegb + 0.33/2. + fgkCellDepth/2.;
  TVirtualMC::GetMC()->Gspos("EHC2", 1, "EUV2", 0., 0.,  zehc2, 0, "ONLY");
  
  //(4)   NEXT PUT THE G10 Box : EDGB
  zedgb = zehc2 + fgkCellDepth/2.+ 0.4/2.;
  TVirtualMC::GetMC()->Gspos("EDGB", 1, "EUV2", 0., 0., zedgb, 0, "ONLY");
  
  
  //===================================================================//
  //---------------------- UM Type B completed ------------------------//
  
}

//_______________________________________________________________________

void AliPMDv1::CreatePMD()
{
  // Create final detector from Unit Modules
  // -- Author : Bedanga and Viyogi June 2003
  
  
  Float_t   zp = fgkZdist;  //Z-distance of PMD from Interaction Point 

  Int_t jhrot12,jhrot13, irotdm;
  Int_t *idtmed = fIdtmed->GetArray()-599;
  
  AliMatrix(irotdm, 90., 0.,  90.,  90., 180., 0.);
  AliMatrix(jhrot12, 90., 180., 90., 270., 0., 0.);
  AliMatrix(jhrot13, 90., 240., 90., 330., 0., 0.);
  
  // Now We Will Calculate Position Co-ordinates of EUM1 & EUV1 in EPM1 & EPM2
  
  Float_t dbox1[3];
  dbox1[0] = fgkCellRadius/fgkSqroot3by2;
  dbox1[1] = fgkNrowUM1*fgkCellRadius;
  dbox1[2] = fgkCellDepth/2.;
  
  Float_t dbox3[3];
  dbox3[0] = (dbox1[0]*fgkNcolUM1)-
    (fgkCellRadius*fgkSqroot3*(fgkNcolUM1-1)/6.) + 0.15 + 0.05 + 0.05;  
  dbox3[1] = dbox1[1]+fgkCellRadius/2. + 0.15 + 0.05 + 0.05; 
  dbox3[2] = fgkCellDepth/2.;
 
  Float_t dboxCGA[3];
  dboxCGA[0]  = dbox3[0] - 0.15; 
  dboxCGA[1]  = dbox3[1] - 0.15; 
  dboxCGA[2]  = 0.4/2.;

  Float_t dboxSS1[3];
  dboxSS1[0]   = dboxCGA[0]+fgkSSBoundary; 
  dboxSS1[1]   = dboxCGA[1]+fgkSSBoundary;       
  dboxSS1[2]   = fgkThSS/2.; 

  Float_t dboxUM1[3];
  dboxUM1[0] = dboxSS1[0];
  dboxUM1[1] = dboxSS1[1];
  dboxUM1[2] = fgkThSS/2. + 1.2;

  Float_t dboxSM1[3];
  dboxSM1[0] = fSMLengthax + 0.05; // 0.05cm for the ESC1,2 
  dboxSM1[1] = fSMLengthay;
  dboxSM1[2] = dboxUM1[2];
 
  // Position co-ordinates of the unit modules in EPM1 & EPM2
  Float_t xa1,xa2,xa3,ya1,ya2; 
  xa1 =  dboxSM1[0] - dboxUM1[0];
  xa2 = xa1 - dboxUM1[0] - 0.1 - dboxUM1[0];
  xa3 = xa2 - dboxUM1[0] - 0.1 - dboxUM1[0];
  ya1 = dboxSM1[1]  - 0.2 - dboxUM1[1];
  ya2 = ya1 - dboxUM1[1] - 0.3 - dboxUM1[1];
  
  // Next to Calculate Position Co-ordinates of EUM2 & EUV2 in EPM3 & EPM4
  
  Float_t dbox2[3];
  dbox2[1] = fgkNrowUM2*fgkCellRadius;
  dbox2[0] = dbox1[0];
  dbox2[2] = dbox1[2];
  
  Float_t dbox4[3];
  dbox4[0] =(dbox2[0]*fgkNcolUM2)-
    (fgkCellRadius*fgkSqroot3*(fgkNcolUM2-1)/6.) + 0.15 + 0.05 + 0.05; 
  dbox4[1] = dbox2[1] + fgkCellRadius/2. + 0.15 + 0.05 + 0.05;
  dbox4[2] = dbox3[2];
  
  Float_t dboxCGB[3];
  dboxCGB[0]  = dbox4[0] - 0.15; 
  dboxCGB[1]  = dbox4[1] - 0.15; 
  dboxCGB[2]  = 0.4/2.;
  
  Float_t dboxSS2[3];
  dboxSS2[0]  = dboxCGB[0] + fgkSSBoundary; 
  dboxSS2[1]  = dboxCGB[1] + fgkSSBoundary;       
  dboxSS2[2]  = fgkThSS/2.;
  
  Float_t dboxUM2[3];
  dboxUM2[0] = dboxSS2[0];
  dboxUM2[1] = dboxSS2[1];
  dboxUM2[2] = fgkThSS/2. + 2.4/2.; // 2.4 cm is added for  FEE Board thickness

  Float_t dboxSM2[3];
  dboxSM2[0] = fSMLengthbx + 0.05;  // 0.05cm for the ESC3,4
  dboxSM2[1] = fSMLengthby;
  dboxSM2[2] = dboxUM2[2];
  
  // Position co-ordinates of the unit modules in EPM3 & EPM4 
  // Space is added to provide a gapping for HV between UM's
  Float_t xb1,xb2,yb1,yb2,yb3; 
  xb1 = dboxSM2[0] - 0.1 - dboxUM2[0];
  xb2 = xb1 - dboxUM2[0] - 0.1 - dboxUM2[0];
  yb1 = dboxSM2[1] -  0.2 - dboxUM2[1];
  yb2 = yb1 - dboxUM2[1] - 0.2 -  dboxUM2[1];
  yb3 = yb2 - dboxUM2[1] - 0.3-  dboxUM2[1];


  // Create Volumes for Lead(Pb) Plates

  // Lead Plate For LONG TYPE
  // X-dimension of Lead Plate = 3*(X-dimension of EUM1 or EUV1) + gap provided between unit modules 
  // Y-dimension of Lead Plate = 2*(Y-dimension of EUM1 or EUV1) + thickness of SS channels 
  // + tolerance
  // Z-demension of Lead Plate = 1.5cm 
  // Integer assigned to Pb-medium is 600

   Float_t dboxLeadA[3];
  dboxLeadA[0] = fSMLengthax; 
  dboxLeadA[1] = fSMLengthay;
  dboxLeadA[2] = fgkThLead/2.;

  TVirtualMC::GetMC()->Gsvolu("ELDA","BOX", idtmed[600], dboxLeadA, 3);

  //LEAD Plate For SHORT TYPE
  // X-dimension of Lead Plate = 2*(X-dimension of EUM2 or EUV2) + gap provided between unit modules 
  // Y-dimension of Lead Plate = 3*(Y-dimension of EUM2 or EUV2) + thickness of SS channels 
  // + tolerance
  // Z-demension of Lead Plate = 1.5cm 
  // Integer assigned to Pb-medium is 600

   Float_t dboxLeadB[3];
  dboxLeadB[0] = fSMLengthbx; 
  dboxLeadB[1] = fSMLengthby; 
  dboxLeadB[2] = fgkThLead/2.;

  TVirtualMC::GetMC()->Gsvolu("ELDB","BOX", idtmed[600], dboxLeadB, 3);

  //=========== CREATE MOTHER VOLUMES FOR PMD ===========================/

  Float_t serviceX    = 23.2;
  Float_t serviceYa   = 5.2;
  Float_t serviceYb   = 9.8;
  Float_t serviceXext = 16.0;

  // Five Mother Volumes of PMD are Created
  // Two Volumes EPM1 & EPM2 of Long Type
  // Other Two Volumes EPM3 & EPM4 for Short Type
  // Fifth Volume EFGD for Girders and its Carriage
  // Four Volmes EPM1, EPM2, EPM3 & EPM4 are Placed such that
  // to create a hole and avoid overlap with Beam Pipe

  // Create Volume FOR EPM1 
  // X-dimension = fSMLengthax + Extended Iron Support(23.2cm) + 
  // Extension in Module(16cm) for full coverage of Detector + 1mm thick SS-Plate
  // Y-dimension = fSMLengthay + Extended Iron Support(5.2cm)
  // Z-dimension = fSMthick/2.; fSMthick=17cm is full profile of PMD in Z-Side
  // Note:- EPM1 is a Volume of Air

  Float_t gaspmd1[3];
  gaspmd1[0] = fSMLengthax + serviceX/2.+ serviceXext/2. + 0.05; //0.05cm for the thickness of 
  gaspmd1[1] = fSMLengthay + serviceYa/2.;                       //SS-plate for cooling encloser  
  gaspmd1[2] = fSMthick/2.;
  
  TVirtualMC::GetMC()->Gsvolu("EPM1", "BOX", idtmed[698], gaspmd1, 3);


  // Create Volume FOR EPM2 

  // X-dimension = fSMLengthax + Extended Iron Support(23.2cm) + 
  // Extension in Module(16cm) for full coverage of Detector + 1mm thick SS-Plate
  // Y-dimension = fSMLengthay + Extended Iron Support(9.8cm)
  // Z-dimension = fSMthick/2.; fSMthick=17cm is full profile of PMD in Z-Side
  // Note:- EPM2 is a Volume of Air

  Float_t gaspmd2[3];
  gaspmd2[0] = fSMLengthax + serviceX/2. + serviceXext/2. + 0.05; //0.05cm for the thickness of 
  gaspmd2[1] = fSMLengthay + serviceYb/2.;                        //SS-plate for cooling encloser
  gaspmd2[2] = fSMthick/2.;

  TVirtualMC::GetMC()->Gsvolu("EPM2", "BOX", idtmed[698], gaspmd2, 3);

  // Create Volume FOR EPM3

  // X-dimension = fSMLengthbx + Extended Iron Support(23.2cm) + 
  // Extension in Module(16cm) for full coverage of Detector
  // Y-dimension = fSMLengthby + Extended Iron Support(5.2cm)
  // Z-dimension = fSMthick/2.; fSMthick=17cm is full profile of PMD in Z-Side
  // Note:- EPM3 is a Volume of Air


  Float_t gaspmd3[3];
  gaspmd3[0] = fSMLengthbx + serviceX/2. + serviceXext/2.+ 0.05; //0.05cm for the thickness of  
  gaspmd3[1] = fSMLengthby + serviceYa/2.;                       //SS-plate for cooling encloser  
  gaspmd3[2] = fSMthick/2.;

  TVirtualMC::GetMC()->Gsvolu("EPM3", "BOX", idtmed[698], gaspmd3, 3);

  // Create Volume FOR EPM4

  // X-dimension = fSMLengthbx + Extended Iron Support(23.2cm) + 
  // Extension in Module(16cm) for full coverage of Detector
  // Y-dimension = fSMLengthby + Extended Iron Support(9.8cm)
  // Z-dimension = fSMthick/2.; fSMthick=17cm is full profile of PMD in Z-Side
  // Note:- EPM4 is a Volume of Air
  
  Float_t gaspmd4[3];
  gaspmd4[0] = fSMLengthbx + serviceX/2. + serviceXext/2.+ 0.05;  //0.05cm for the thickness of
  gaspmd4[1] = fSMLengthby + serviceYb/2.;                        //SS-plate for cooling encloser   
  gaspmd4[2] = fSMthick/2.;

  TVirtualMC::GetMC()->Gsvolu("EPM4", "BOX", idtmed[698], gaspmd4, 3);
  
  //  Create the Fifth Mother Volume of Girders and its Carriage
  //-------------------------------------------------------------//
  // Create the Girders
  
  // X-dimension = 238.7cm 
  // Y-dimension = 12.0cm 
  // Z-dimension = 7.0cm 
  // Girders are the Volume of Iron
  // And the Integer Assigned to SS is 618

  Float_t grdr[3];
  grdr[0] = 238.7/2.;
  grdr[1] = 12.0/2.;
  grdr[2] = 7.0/2.; 

  TVirtualMC::GetMC()->Gsvolu("EGDR", "BOX", idtmed[618], grdr, 3);
 
  // Create Air Strip for Girders as the Girders are hollow
  // Girders are 1cm thick in Y and Z on both sides
 
  Float_t airgrdr[3];
  airgrdr[0] = grdr[0];
  airgrdr[1] = grdr[1] - 1.0;
  airgrdr[2] = grdr[2] - 1.0;
  
  TVirtualMC::GetMC()->Gsvolu("EAIR", "BOX", idtmed[698], airgrdr, 3);

  // Positioning the air strip EAIR in girder EGDR  
  TVirtualMC::GetMC()->Gspos("EAIR", 1, "EGDR",  0., 0., 0.,  0, "ONLY");
  
  // Create the Carriage for Girders
  // Originally, Carriage is divided in two parts
  // 64.6cm on -X side, 44.2cm on +X side and 8.2cm is the gap between two
  // In approximation we have taken these together as a single Volume
  // With X = 64.6cm + 44.2cm + 8.2cm
  // Y-dimension = 4.7cm
  // Z-dimension = 18.5cm
  // Carriage is a Volume of SS
    
  Float_t xgrdr[3];
  xgrdr[0] = (64.6 + 44.2 + 8.2)/2.;  
  xgrdr[1] = 4.7/2.; 
  xgrdr[2] = 18.5/2.;

  TVirtualMC::GetMC()->Gsvolu("EXGD", "BOX", idtmed[618], xgrdr, 3);

  // Create Air Strip for the Carriage EXGD as it is hollow
  // Carriage is 1cm thick in Y on one side and in Z on both sides 

  Float_t xairgrdr[3];
  xairgrdr[0] = xgrdr[0];
  xairgrdr[1] = xgrdr[1] - 0.5;
  xairgrdr[2] = xgrdr[2] - 1.0;
  
  TVirtualMC::GetMC()->Gsvolu("EXIR", "BOX", idtmed[698], xairgrdr, 3);
  
  // Positioning the air strip EXIR in CArriage EXGD
  TVirtualMC::GetMC()->Gspos("EXIR", 1, "EXGD",  0., -0.05, 0.,  0, "ONLY");

  // Now Create the master volume of air containing Girders & Carriage
    
  // X-dimension = same as X-dimension of Girders(EGDR)
  // Y-dimension = Y of Girder(EGDR) + Y of Carriage(EXGD) + gap between two
  // Z-dimenson = same as Z of Carriage(EXGD)
  // Note:- It is a volume of Air

  Float_t fulgrdr[3];
  fulgrdr[0] = 238.7/2.;
  fulgrdr[1] = 17.5/2.; 
  fulgrdr[2] = 18.5/2.;

  TVirtualMC::GetMC()->Gsvolu("EFGD", "BOX", idtmed[698], fulgrdr, 3);

  // Positioning the EGDR and EXGD in EFGD

  TVirtualMC::GetMC()->Gspos("EXGD", 1, "EFGD",  0., 6.4, 0.,      0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EGDR", 1, "EFGD",  0., -2.75, -5.75, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EGDR", 2, "EFGD",  0., -2.75, 5.75,  0, "ONLY");

  //=========== Mother Volumes are Created ============================//

  // Create the Volume of 1mm thick SS-Plate  for cooling encloser
  // These are placed on the side close to the Beam Pipe
  // SS-Plate is perpendicular to the plane of Detector 
 
  // For LONG TYPE

  // For EPM1
  // X-dimension = 0.1cm
  // Y-dimension = same as Y of EPM1
  // Z-dimension = Y of EPM1 - 0.1; 0.1cm is subtracted as 1mm thick 
  // FR4 sheets for the detector encloser placed on both sides
  // It is a Volume of SS
  // Integer assigned to SS is 618
 
  Float_t sscoolencl1[3];
  sscoolencl1[0] = 0.05;  
  sscoolencl1[1] = gaspmd1[1];
  sscoolencl1[2] = gaspmd1[2] - 0.2/2.;

  TVirtualMC::GetMC()->Gsvolu("ESC1", "BOX", idtmed[618], sscoolencl1, 3);

  // Placement of ESC1  in EPM1
  TVirtualMC::GetMC()->Gspos("ESC1", 1,  "EPM1", -gaspmd1[0] + 0.05, 0., 0., 0, "ONLY");


  // For EPM2
  // X-dimension = 0.1cm
  // Y-dimension = same as Y of EPM2
  // Z-dimension = Y of EPM2 - 0.1; 0.1cm is subtracted as 1mm thick 
  // FR4 sheets for the detector encloser placed on both sides
  // It is a Volume of SS
 
  Float_t sscoolencl2[3];
  sscoolencl2[0] = 0.05;  
  sscoolencl2[1] = gaspmd2[1];
  sscoolencl2[2] = gaspmd2[2] - 0.2/2.;

  TVirtualMC::GetMC()->Gsvolu("ESC2", "BOX", idtmed[618], sscoolencl2, 3);

  // Placement of ESC2  in EPM2
  TVirtualMC::GetMC()->Gspos("ESC2", 1,  "EPM2",    gaspmd2[0] - 0.05 , 0., 0., 0, "ONLY");

  // For SHORT TYPE

  // For EPM3
  // X-dimension = 0.1cm
  // Y-dimension = same as Y of EPM3
  // Z-dimension = Y of EPM3 - 0.1; 0.1cm is subtracted as 1mm thick 
  // FR4 sheets for the detector encloser placed on both sides
  // It is a Volume of SS
  
  Float_t sscoolencl3[3];
  sscoolencl3[0] = 0.05;  
  sscoolencl3[1] = gaspmd3[1];
  sscoolencl3[2] = gaspmd3[2] - 0.2/2.;

  TVirtualMC::GetMC()->Gsvolu("ESC3", "BOX", idtmed[618], sscoolencl3, 3);

  // Placement of ESC3  in EPM3
  TVirtualMC::GetMC()->Gspos("ESC3", 1,  "EPM3",    gaspmd3[0] - 0.05 , 0., 0., 0, "ONLY");


  // For EPM4
  // X-dimension = 0.1cm
  // Y-dimension = same as Y of EPM4
  // Z-dimension = Y of EPM4 - 0.1; 0.1cm is subtracted as 1mm thick 
  // FR4 sheets for the detector encloser placed on both sides
  // It is a Volume of SS
 
  Float_t sscoolencl4[3];
  sscoolencl4[0] = 0.05;  
  sscoolencl4[1] = gaspmd4[1];
  sscoolencl4[2] = gaspmd4[2] - 0.2/2.;

  TVirtualMC::GetMC()->Gsvolu("ESC4", "BOX", idtmed[618], sscoolencl4, 3);

  // Placement of ESC4  in EPM4
  TVirtualMC::GetMC()->Gspos("ESC4", 1, "EPM4", -gaspmd4[0] + 0.05 , 0., 0., 0, "ONLY");

  //======== CREATE SS SUPPORTS FOR EPM1, EPM2, EPM3 & EPM4 =========//
  // --- DEFINE SS volumes  for EPM1 & EPM2 ---

  // Create SS Support For EPM1

  // X-dimension = fSMLengthax + Extended Iron Support(23.2cm)
  // Y-dimension = fSMLengthay + Extended Iron Support(5.2cm)
  // Z-dimension = thickness of Iron support(0.5cm)
  // It is a Volume of SS
  // Integer assigned to SS is 618

  Float_t dboxFea1[3];
  dboxFea1[0] = fSMLengthax + serviceX/2.;  
  dboxFea1[1] = fSMLengthay + serviceYa/2.;
  dboxFea1[2] = fgkThSteel/2.;
  
  TVirtualMC::GetMC()->Gsvolu("EFE1","BOX", idtmed[618], dboxFea1, 3);


  // Create SS Support For EPM2

  // X-dimension = fSMLengthax + Extended Iron Support(23.2cm)
  // Y-dimension = fSMLengthay + Extended Iron Support(9.8cm)
  // Z-dimension = thickness of Iron support(0.5cm)
  // It is a Volume of SS
  // Integer assigned to SS is 618

  Float_t dboxFea2[3];
  dboxFea2[0] = fSMLengthax + serviceX/2.;   
  dboxFea2[1] = fSMLengthay + serviceYb/2.;  
  dboxFea2[2] = fgkThSteel/2.;
  
  TVirtualMC::GetMC()->Gsvolu("EFE2","BOX", idtmed[618], dboxFea2, 3);

  // Create SS Support For EPM3

  // X-dimension = fSMLengthbx + Extended Iron Support(23.2cm)
  // Y-dimension = fSMLengthby + Extended Iron Support(5.2cm)
  // Z-dimension = thickness of Iron support(0.5cm)
  // It is a Volume of SS
  // Integer assigned to SS is 618  

  Float_t dboxFea3[3];
  dboxFea3[0] = fSMLengthbx + serviceX/2.; 
  dboxFea3[1] = fSMLengthby + serviceYa/2.;
  dboxFea3[2] = fgkThSteel/2.;
  
  TVirtualMC::GetMC()->Gsvolu("EFE3","BOX", idtmed[618], dboxFea3, 3);

  // Create SS Support For EPM4

  // X-dimension = fSMLengthbx + Extended Iron Support(23.2cm)
  // Y-dimension = fSMLengthby + Extended Iron Support(9.8cm)
  // Z-dimension = thickness of Iron support(0.5cm)
  // It is a Volume of SS
  // Integer assigned to SS is 618  
 
  Float_t dboxFea4[3];
  dboxFea4[0] = fSMLengthbx + serviceX/2.;  
  dboxFea4[1] = fSMLengthby + serviceYb/2.; 
  dboxFea4[2] = fgkThSteel/2.;
  
  TVirtualMC::GetMC()->Gsvolu("EFE4","BOX", idtmed[618], dboxFea4, 3);


  //=============== Volumes for SS support are Completed =============//

  // Create FR4 Sheets to enclose the PMD which are Placed parallel to the
  // plane of the detector. Four FR4 sheets are created with the dimensions
  // corresponding to the Iron Supports
  // This is cooling encloser.

  // Create FR4 sheet ECC1
  // X-dimension = same as EFE1
  // Y-dimension = same as EFE1
  // Z-dimension = 0.1cm
  // FR4 medium is same as that of G10
  // Integer assigned to FR4 medium is 607

  Float_t enclos1[3];
  enclos1[0] = dboxFea1[0];   
  enclos1[1] = dboxFea1[1];
  enclos1[2] = 0.05;

  TVirtualMC::GetMC()->Gsvolu("ECC1", "BOX", idtmed[607], enclos1, 3);

  // Create FR4 sheet ECC2
  // X-dimension = same as EFE2
  // Y-dimension = same as EFE2
  // Z-dimension = 0.1cm

  Float_t enclos2[3];
  enclos2[0] = dboxFea2[0];  
  enclos2[1] = dboxFea2[1];
  enclos2[2] = 0.05;

  TVirtualMC::GetMC()->Gsvolu("ECC2", "BOX", idtmed[607], enclos2, 3);

  // Create FR4 sheet ECC3
  // X-dimension = same as EFE3
  // Y-dimension = same as EFE3
  // Z-dimension = 0.1cm

  Float_t enclos3[3];
  enclos3[0] = dboxFea3[0];  
  enclos3[1] = dboxFea3[1];
  enclos3[2] = 0.05;
  
  TVirtualMC::GetMC()->Gsvolu("ECC3", "BOX", idtmed[607], enclos3, 3);
  
  // Create FR4 sheet ECC4
  // X-dimension = same as EFE4
  // Y-dimension = same as EFE4
  // Z-dimension = 0.1cm

  Float_t enclos4[3];
  enclos4[0] = dboxFea4[0];   
  enclos4[1] = dboxFea4[1];
  enclos4[2] = 0.05;

  TVirtualMC::GetMC()->Gsvolu("ECC4", "BOX", idtmed[607], enclos4, 3);

  //--------------- FR4 SHEETS COMPLETED ---------------------------//

  //------------- Create the SS-Channels(Horizontal Rails) to Place
  //     Unit Modules on SS Support -------------------------------------//
  
  // Two types of SS-Channels are created 
  // as we have two types of modules
  
  // Create SS-channel for Long Type
  // X-dimension = same as Lead Plate ELDA
  // Y-dimension = 0.1cm
  // Z-dimension = 2.0cm
  // Volume medium is SS

  Float_t channel12[3];
  channel12[0] = fSMLengthax;  
  channel12[1] = 0.05; 
  channel12[2] = 2.0/2.; 

  TVirtualMC::GetMC()->Gsvolu("ECHA", "BOX", idtmed[618], channel12, 3);
  
  // Create SS-channel for Short Type
  // X-dimension = same as Lead Plate ELDB
  // Y-dimension = 0.1cm
  // Z-dimension = 2.0cm
  // Volume medium is SS

  Float_t channel34[3];
  channel34[0] = fSMLengthbx;  
  channel34[1] = 0.05; 
  channel34[2] = 2.0/2.; 

  TVirtualMC::GetMC()->Gsvolu("ECHB", "BOX", idtmed[618], channel34, 3);

  //----------------- SS-Channels are Copmleted --------------------//

  //========= POSITIONING OF SS SUPPORT AND LEAD PLATES IN QUADRANTS =====//
  
  /**************** Z-Distances of different Components **********/
  
  Float_t zcva,zfea,zpba,zpsa,zchanVeto,zchanPS, zelvdbVeto, zelvdbPS;
  
  
  zpba       =  - fgkThSteel/2.;                         //z-position of Pb plate
  zfea       =  fgkThLead/2.;                            //z-position of SS-Support
  zchanVeto  =  zpba -  fgkThLead/2. - channel12[2];     //z-position of SS-channel on Veto
  zchanPS    =  zfea + fgkThSteel/2. + channel12[2];     //z-position of SS-channel on Preshower
  zpsa       =  zfea + fgkThSteel/2. + fDthick;          //z-position of Preshower
  zcva       =  zpba - fgkThLead/2.- fDthick;            //z-position of Veto
  
  zelvdbVeto =  zpba + fgkThLead/2.  - 8.9/2.;           //z-position of LVDBs on Veto side
  zelvdbPS   =  zfea + fgkThSteel/2. + 7.4/2.;           //z-position of LVDBs on Preshower side
  
  // FOR LONG TYPE
  Float_t  xLead1,yLead1,zLead1, xLead2,yLead2,zLead2;
  Float_t  xIron1,yIron1,zIron1, xIron2,yIron2,zIron2;
  
  
  xIron1 = - 16.0/2. + 0.1/2.; // half of 0.1cm is added as 1mm SS sheet is placed 
  yIron1 = 0.;
  zIron1 = zfea;
  
  xIron2 = 16.0/2. - 0.1/2.;  // half of 0.1cm is added as 1mm SS sheet is placed 
  yIron2 = 0.;
  zIron2 = zfea;    

  
  xLead1 = xIron1 - 23.2/2.; 
  yLead1 = -5.2/2.;
  zLead1 = zpba;
  
  xLead2 =xIron2 + 23.2/2.; 
  yLead2 = 9.8/2.;
  zLead2 = zpba;    
  
  TVirtualMC::GetMC()->Gspos("EFE1", 1, "EPM1", xIron1,  yIron1, zfea, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ELDA", 1, "EPM1", xLead1,  yLead1, zpba, 0, "ONLY"); 
  TVirtualMC::GetMC()->Gspos("EFE2", 1, "EPM2", xIron2,  yIron2, zfea, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ELDA", 1, "EPM2", xLead2,  yLead2, zpba, jhrot12, "ONLY"); 
  
  
  // FOR SHORT TYPE
  Float_t xLead3,yLead3,zLead3, xLead4,yLead4,zLead4;
  Float_t xIron3,yIron3,zIron3, xIron4,yIron4,zIron4;
  
  
  xIron3 =  16.0/2.- 0.1/2.;  // half of 0.1cm is added as 1mm SS sheet is placed ; 
  yIron3 = 0.;
  zIron3 = zfea;
  
  xIron4 = - 16.0/2.+ 0.1/2.; // half of 0.1cm is added as 1mm SS sheet is placed; 
  yIron4 = 0.;
  zIron4 = zfea;    
  
  xLead3 = xIron3 + 23.2/2.; 
  yLead3 = -5.2/2.;
  zLead3 = zpba;
  
  xLead4 = xIron4 - 23.2/2.; 
  yLead4 = 9.8/2.;
  zLead4 = zpba;    
  
  TVirtualMC::GetMC()->Gspos("EFE3", 1,  "EPM3",  xIron3,  yIron3,  zfea, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ELDB", 1,  "EPM3",  xLead3,  yLead3,  zpba, 0, "ONLY"); 
  TVirtualMC::GetMC()->Gspos("EFE4", 1,  "EPM4",  xIron4,  yIron4,  zfea, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ELDB", 1,  "EPM4",  xLead4,  yLead4,  zpba, jhrot12, "ONLY"); 
  
  //===================================================================//
  // Placement of FR4 sheets as encloser of full profile of PMD

  TVirtualMC::GetMC()->Gspos("ECC1", 1, "EPM1",  xIron1, yIron1, -8.45,  0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECC2", 1, "EPM2",  xIron2, yIron2, -8.45,  0,"ONLY");
  TVirtualMC::GetMC()->Gspos("ECC3", 1, "EPM3",  xIron3, yIron3, -8.45, 0,"ONLY");
  TVirtualMC::GetMC()->Gspos("ECC4", 1, "EPM4",  xIron4, yIron4, -8.45, 0,"ONLY");

  TVirtualMC::GetMC()->Gspos("ECC1", 2, "EPM1",  xIron1, yIron1,  8.45, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECC2", 2, "EPM2",  xIron2, yIron2,  8.45, 0,"ONLY");
  TVirtualMC::GetMC()->Gspos("ECC3", 2, "EPM3",  xIron3, yIron3,  8.45, 0,"ONLY");
  TVirtualMC::GetMC()->Gspos("ECC4", 2, "EPM4",  xIron4, yIron4,  8.45, 0,"ONLY");

  //----------------- NOW TO PLACE SS-CHANNELS -----------------------// 
  
  Float_t xchanepm11, ychanepm11,ychanepm12;
  Float_t xchanepm21, ychanepm21,ychanepm22;
  Float_t xchanepm31, ychanepm31,ychanepm32,ychanepm33,ychanepm34;
  Float_t xchanepm41, ychanepm41,ychanepm42,ychanepm43,ychanepm44;
  
  xchanepm11 = xLead1;
  ychanepm11 = ya1 + yLead1 + dboxSS1[1] + 0.1 + 0.1/2.;
  ychanepm12 = ya1 + yLead1 - dboxSS1[1] - 0.1 - 0.1/2.;
  
  xchanepm21 = xLead2;
  ychanepm21 = -ya1 + yLead2 - dboxSS1[1] - 0.1 - 0.1/2.;
  ychanepm22 = -ya1 + yLead2 + dboxSS1[1] + 0.1 + 0.1/2.;
  
  TVirtualMC::GetMC()->Gspos("ECHA", 1, "EPM1", xchanepm11, ychanepm11, zchanPS,   0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECHA", 2, "EPM1", xchanepm11, ychanepm12, zchanPS,   0, "ONLY"); 
  TVirtualMC::GetMC()->Gspos("ECHA", 3, "EPM1", xchanepm11, ychanepm11, zchanVeto, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECHA", 4, "EPM1", xchanepm11, ychanepm12, zchanVeto, 0, "ONLY"); 
  TVirtualMC::GetMC()->Gspos("ECHA", 1, "EPM2", xchanepm21, ychanepm21, zchanPS,   0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECHA", 2, "EPM2", xchanepm21, ychanepm22, zchanPS,   0, "ONLY"); 
  TVirtualMC::GetMC()->Gspos("ECHA", 3, "EPM2", xchanepm21, ychanepm21, zchanVeto, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECHA", 4, "EPM2", xchanepm21, ychanepm22, zchanVeto, 0, "ONLY"); 
  
  xchanepm31 = xLead3;
  ychanepm31 = yb1 + yLead3 + dboxSS2[1] + 0.1 + 0.1/2.;
  ychanepm32 = yb1 + yLead3 - dboxSS2[1] - 0.1 - 0.1/2.;
  ychanepm33 = yb3 + yLead3 + dboxSS2[1] + 0.1 + 0.1/2.;
  ychanepm34 = yb3 + yLead3 - dboxSS2[1] - 0.1 - 0.1/2.;
  
  xchanepm41 = xLead4;
  ychanepm41 = -yb1 + yLead4 - dboxSS2[1] - 0.1 - 0.1/2.;
  ychanepm42 = -yb1 + yLead4 + dboxSS2[1] + 0.1 + 0.1/2.;
  ychanepm43 = -yb3 + yLead4 - dboxSS2[1] - 0.1 - 0.1/2.;
  ychanepm44 = -yb3 + yLead4 + dboxSS2[1] + 0.1 + 0.1/2.;
  
  
  TVirtualMC::GetMC()->Gspos("ECHB", 1, "EPM3", xchanepm31, ychanepm31, zchanPS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECHB", 2, "EPM3", xchanepm31, ychanepm32, zchanPS, 0, "ONLY"); 
  TVirtualMC::GetMC()->Gspos("ECHB", 3, "EPM3", xchanepm31, ychanepm33, zchanPS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECHB", 4, "EPM3", xchanepm31, ychanepm34 + 0.200005, zchanPS, 0, "ONLY"); 
  // Because of overlaping a factor 0.200005 is added in ychanepm34
  
  TVirtualMC::GetMC()->Gspos("ECHB", 5, "EPM3", xchanepm31, ychanepm31, zchanVeto, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECHB", 6, "EPM3", xchanepm31, ychanepm32, zchanVeto, 0, "ONLY"); 
  TVirtualMC::GetMC()->Gspos("ECHB", 7, "EPM3", xchanepm31, ychanepm33, zchanVeto, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECHB", 8, "EPM3", xchanepm31, ychanepm34 + 0.200005, zchanVeto, 0, "ONLY"); 
  // Because of overlaping a factor 0.200005 is added in ychanepm34
  
  TVirtualMC::GetMC()->Gspos("ECHB", 1, "EPM4", xchanepm41, ychanepm41, zchanPS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECHB", 2, "EPM4", xchanepm41, ychanepm42, zchanPS, 0, "ONLY"); 
  TVirtualMC::GetMC()->Gspos("ECHB", 3, "EPM4", xchanepm41, ychanepm43, zchanPS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECHB", 4, "EPM4", xchanepm41, ychanepm44 - 0.200002, zchanPS, 0, "ONLY"); 
  // Because of overlaping a factor 0.200002 is subtracted in ychanepm44

  TVirtualMC::GetMC()->Gspos("ECHB", 5, "EPM4", xchanepm41, ychanepm41, zchanVeto, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECHB", 6, "EPM4", xchanepm41, ychanepm42, zchanVeto, 0, "ONLY"); 
  TVirtualMC::GetMC()->Gspos("ECHB", 7, "EPM4", xchanepm41, ychanepm43, zchanVeto, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECHB", 8, "EPM4", xchanepm41, ychanepm44 -0.200002, zchanVeto, 0, "ONLY"); 
  // Because of overlaping a factor 0.200002 is subtracted in ychanepm44

  //================= Channel Placement Completed  ======================//
  //============ Now to Create Al Box and then LVDBs and Cables          //
  //             are Placed inside it                                    //

  // Eight Al Boxes are created, four on Preshower side 
  // and four on Veto side

  // FOR PRESHOWER

  // First to Create hollow Al Box
  // there are two types of modules, therefore, two Al box of
  // long type and two of short type are created

  // For Long Type
  // X-dimension = 16.5cm
  // Y-dimension = same as EFE1
  // Z-dimension = 7.4cm
  // Integer assigned to Al medium is 603

  Float_t esvdA1[3];
  esvdA1[0]= 16.5/2.;
  esvdA1[1]= dboxFea1[1];
  esvdA1[2]= 7.4/2.;
  
  TVirtualMC::GetMC()->Gsvolu("ESV1", "BOX", idtmed[603], esvdA1, 3);
  TVirtualMC::GetMC()->Gsvolu("ESV2", "BOX", idtmed[603], esvdA1, 3);
  
  // Create Air strip for Al Boxes type-A
  // Al boxes are 3mm thick In X and Z on both sides
  // X-dimension = 16.5cm - 0.3cm
  // Y-dimension = same as EFE1
  // Z-dimension = 7.4cm - 0.3cm

  Float_t eairA1[3];
  eairA1[0]= esvdA1[0] - 0.3;
  eairA1[1]= esvdA1[1];
  eairA1[2]= esvdA1[2] - 0.3;

  TVirtualMC::GetMC()->Gsvolu("EIR1", "BOX", idtmed[698], eairA1, 3);
  TVirtualMC::GetMC()->Gsvolu("EIR2", "BOX", idtmed[698], eairA1, 3);

  // Put air strips EIR1 & EIR2 inside ESV1 & ESV2 respectively    
  TVirtualMC::GetMC()->Gspos("EIR1", 1,  "ESV1", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EIR2", 1,  "ESV2", 0., 0., 0., 0, "ONLY");
  

  // For Short Type
  // X-dimension = 16.5cm
  // Y-dimension = same as EFE3
  // Z-dimension = 7.4cm
  
  Float_t esvdA2[3];
  esvdA2[0]= esvdA1[0];
  esvdA2[1]= dboxFea3[1];
  esvdA2[2]= esvdA1[2];

  TVirtualMC::GetMC()->Gsvolu("ESV3", "BOX", idtmed[603], esvdA2, 3);
  TVirtualMC::GetMC()->Gsvolu("ESV4", "BOX", idtmed[603], esvdA2, 3);
  
  // Create Air strip for Al Boxes type-B
  // Al boxes are 3mm thick In X and Z on both sides
  // X-dimension = 16.5cm - 0.3cm
  // Y-dimension = same as EFE3
  // Z-dimension = 7.4cm - 0.3cm

  Float_t eairA2[3];
  eairA2[0]= esvdA2[0] - 0.3;
  eairA2[1]= esvdA2[1];
  eairA2[2]= esvdA2[2] - 0.3;

  TVirtualMC::GetMC()->Gsvolu("EIR3", "BOX", idtmed[698], eairA2, 3);
  TVirtualMC::GetMC()->Gsvolu("EIR4", "BOX", idtmed[698], eairA2, 3);
  
  // Put air strips EIR3 & EIR4 inside ESV3 & ESV4 respectively        
  TVirtualMC::GetMC()->Gspos("EIR3", 1,  "ESV3", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EIR4", 1,  "ESV4", 0., 0., 0., 0, "ONLY");
  
  
  // FOR VETO

  // First to Create hollow Al Box
  // there are two types of modules, therefore, two Al box of
  // long type and two of short type are created

  // For Long Type
  // X-dimension = 16.5cm
  // Y-dimension = same as EFE1
  // Z-dimension = 8.9cm
  // Integer assigned to Al medium is 603
  
  Float_t esvdB1[3];
  esvdB1[0]= 16.5/2.;
  esvdB1[1]= dboxFea1[1];
  esvdB1[2]= 8.9/2.;

  TVirtualMC::GetMC()->Gsvolu("EVV1", "BOX", idtmed[603], esvdB1, 3);
  TVirtualMC::GetMC()->Gsvolu("EVV2", "BOX", idtmed[603], esvdB1, 3);

  // Create Air strip for Al Boxes long type
  // Al boxes are 3mm thick In X and Z on both sides
  // X-dimension = 16.5cm - 0.3cm
  // Y-dimension = same as EFE1
  // Z-dimension = 8.9cm - 0.3cm

  Float_t eairB1[3];
  eairB1[0]= esvdB1[0] - 0.3;
  eairB1[1]= esvdB1[1];
  eairB1[2]= esvdB1[2] - 0.3;

  TVirtualMC::GetMC()->Gsvolu("EIR5", "BOX", idtmed[698], eairB1, 3);
  TVirtualMC::GetMC()->Gsvolu("EIR6", "BOX", idtmed[698], eairB1, 3);
 
  // Put air strips EIR5 & EIR6 inside EVV1 & EVV2 respectively        
  TVirtualMC::GetMC()->Gspos("EIR5", 1,  "EVV1", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EIR6", 1,  "EVV2", 0., 0., 0., 0, "ONLY");


  // For Short Type
  // X-dimension = 16.5cm
  // Y-dimension = same as EFE3
  // Z-dimension = 8.9cm
  // Integer assigned to Al medium is 603
  
  Float_t esvdB2[3];
  esvdB2[0]= esvdB1[0];
  esvdB2[1]= dboxFea3[1];
  esvdB2[2]= esvdB1[2];

  TVirtualMC::GetMC()->Gsvolu("EVV3", "BOX", idtmed[603], esvdB2, 3);
  TVirtualMC::GetMC()->Gsvolu("EVV4", "BOX", idtmed[603], esvdB2, 3);

  
  // Create Air strip for Al Boxes short type
  // Al boxes are 3mm thick In X and Z on both sides
  // X-dimension = 16.5cm - 0.3cm
  // Y-dimension = same as EFE3
  // Z-dimension = 8.9cm - 0.3cm
  
  Float_t eairB2[3];
  eairB2[0]= esvdB2[0] - 0.3;
  eairB2[1]= esvdB2[1];
  eairB2[2]= esvdB2[2] - 0.3;
  
  TVirtualMC::GetMC()->Gsvolu("EIR7", "BOX", idtmed[698], eairB2, 3);
  TVirtualMC::GetMC()->Gsvolu("EIR8", "BOX", idtmed[698], eairB2, 3);
  
  // Put air strips EIR7 & EIR8 inside EVV3 & EVV4 respectively      
  TVirtualMC::GetMC()->Gspos("EIR7", 1,  "EVV3", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EIR8", 1,  "EVV4", 0., 0., 0., 0, "ONLY");
  
  //------------ Al Boxes Completed ----------------------/
  
  //--------------Now Create LVDBs----------------------/
  
  // LVDBs are the volumes of G10
  // X-dimension = 10.0cm
  // Y-dimension = 8.0cm
  // Z-dimension = 0.2cm
  // Integer assigned to the G10 medium is 607
  
  Float_t elvdb[3];
  elvdb[0]= 10.0/2.;
  elvdb[1]= 8.0/2.;
  elvdb[2]= 0.2/2.;
  
  TVirtualMC::GetMC()->Gsvolu("ELVD", "BOX", idtmed[607], elvdb, 3);
  

  // Put the LVDBs inside Air Boxes
  Float_t yesvd = dboxFea1[1] - 25.0 - 4.0;
  
  for(Int_t jj =1; jj<=6; jj++){
    
    TVirtualMC::GetMC()->Gspos("ELVD", jj,  "EIR1", 0., yesvd, 0., 0, "ONLY");
    TVirtualMC::GetMC()->Gspos("ELVD", jj,  "EIR2", 0., yesvd, 0., 0, "ONLY");

    yesvd = yesvd -  4.0 - 0.5 - 4.0;
    
  }
  
  yesvd = dboxFea3[1] - 15.0 - 4.0;
  
  for(Int_t jj =1; jj<=6; jj++){
    
    TVirtualMC::GetMC()->Gspos("ELVD", jj,  "EIR3", 0., yesvd, 0., 0, "ONLY");
    TVirtualMC::GetMC()->Gspos("ELVD", jj,  "EIR4", 0., yesvd, 0., 0, "ONLY");

    yesvd = yesvd -  4.0 - 0.5 - 4.0;
  }
  
  yesvd = dboxFea1[1] - 25.0 - 4.0;
  
  for(Int_t jj =1; jj<=6; jj++){
    
    TVirtualMC::GetMC()->Gspos("ELVD", jj,  "EIR5", 0., yesvd, 0., 0, "ONLY");
    TVirtualMC::GetMC()->Gspos("ELVD", jj,  "EIR6", 0., yesvd, 0., 0, "ONLY");

    yesvd = yesvd -  4.0 - 0.5 - 4.0;
  }
  
  yesvd = dboxFea3[1] - 15.0 - 4.0;
  
  for(Int_t jj =1; jj<=6; jj++){
    
    TVirtualMC::GetMC()->Gspos("ELVD", jj,  "EIR7", 0., yesvd, 0., 0, "ONLY");
    TVirtualMC::GetMC()->Gspos("ELVD", jj,  "EIR8", 0., yesvd, 0., 0, "ONLY");

    yesvd = yesvd -  4.0 - 0.5 - 4.0;
  }

  
  //----------------- LVDBs Placement Completed--------------//
  
  // ------------ Now Create Cables ------------------------//
  
  // There are a number of cables
  // We have reduced the number of volumes to 4
  // And these 4 Volumes of Cables are placed repeatedly
  // in the four quadrants (EPM1,2,3,4)
  // The placement of Cables are in good approximations 
  // The material medium for Cables is a mixture of Plastic
  // and Copper(Cu). Therefore, in a good approximation a mixture
  // is created and Integer assigned to this medium is 631
  
  Float_t cable1[3];
  cable1[0] = 2.5/2.;
  cable1[1] = dboxFea1[1];
  cable1[2] = 2.4/2.;
  
  TVirtualMC::GetMC()->Gsvolu("ECB1", "BOX", idtmed[631], cable1, 3);
  
  Float_t cable2[3];
  cable2[0] = 2.5/2.;
  cable2[1] = dboxFea3[1];
  cable2[2] = 2.4/2.;
  
  TVirtualMC::GetMC()->Gsvolu("ECB2", "BOX", idtmed[631], cable2, 3);
  
  Float_t cable3[3];
  cable3[0] = 2.5/2.;
  cable3[1] = dboxFea3[1] - dboxUM2[1];
  cable3[2] = 2.4/2.;
  
  TVirtualMC::GetMC()->Gsvolu("ECB3", "BOX", idtmed[631], cable3, 3);
  
  Float_t cable4[3];
  cable4[0] = 2.5/2.;
  cable4[1] = dboxUM2[1];
  cable4[2] = 2.4/2.;
  
  TVirtualMC::GetMC()->Gsvolu("ECB4", "BOX", idtmed[631], cable4, 3);
  
  // Calculation of the co-ordinates of Cables

  Float_t xcable11pm2, xcable12pm2, xcable2pm1, xcable2pm2,  xcable21pm4,  xcable22pm4;
  Float_t xcable3pm1, xcable3pm3, xcable3pm4, xcable4pm3;

  Float_t ycable2pm1, ycable2pm2;
  Float_t ycable3pm1, ycable3pm3, ycable3pm4, ycable4pm3;
  
  Float_t zcablePS, zcableVeto;
  
  xcable2pm1 = esvdA1[0] - 3.0 - cable1[0];
  xcable3pm1 = xcable2pm1 - cable1[0] - 0.5 -  cable1[0];
  
  xcable11pm2 = -esvdA1[0]+ 3.0 + cable1[0];
  xcable12pm2 = xcable11pm2 + cable1[0] + 0.5 + cable1[0];
  xcable2pm2  = xcable12pm2 + cable1[0] + 0.5 + cable1[0];
  
  xcable3pm3 = -esvdB1[0] + 3.0 + cable1[0];
  xcable4pm3 = xcable3pm3 + cable1[0] + 0.5 + cable1[0];
  
  xcable21pm4 = esvdB1[0] - 3.0 - cable1[0];
  xcable22pm4 = xcable21pm4 - cable1[0] -0.5 - cable1[0];
  xcable3pm4  = xcable22pm4 - cable1[0] -0.5 -cable1[0];
  
  ycable2pm1 = -(esvdA1[1] - esvdA2[1]);
  ycable3pm1 = -esvdA1[1] + cable3[1];
  
  ycable2pm2 =  -(esvdA1[1] - esvdA2[1]);
  
  ycable3pm3 = -dboxUM2[1];
  ycable4pm3 = -esvdA2[1] + dboxUM2[1];
  
  ycable3pm4 = -dboxUM2[1];
  
  zcablePS   = -esvdA1[2] + 0.3 + cable1[2];
  zcableVeto =  esvdB1[2] - 0.3 - cable1[2];



  // Placement of Cables in Air Boxes
  TVirtualMC::GetMC()->Gspos("ECB2", 1,  "EIR1", xcable2pm1, ycable2pm1, zcablePS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB3", 1,  "EIR1", xcable3pm1, ycable3pm1, zcablePS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB2", 1,  "EIR5", xcable2pm1, ycable2pm1, zcableVeto, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB3", 1,  "EIR5", xcable3pm1, ycable3pm1, zcableVeto, 0, "ONLY");
  
  TVirtualMC::GetMC()->Gspos("ECB1", 1,  "EIR2", xcable11pm2,    0.,     zcablePS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB1", 2,  "EIR2", xcable12pm2,    0.,     zcablePS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB2", 1,  "EIR2", xcable2pm2, ycable2pm2, zcablePS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB1", 1,  "EIR6", xcable11pm2,    0.,     zcableVeto, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB1", 2,  "EIR6", xcable12pm2,    0.,     zcableVeto, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB2", 1,  "EIR6", xcable2pm2, ycable2pm2, zcableVeto, 0, "ONLY");
  
  TVirtualMC::GetMC()->Gspos("ECB3", 1,  "EIR3", xcable3pm3, ycable3pm3, zcablePS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB4", 1,  "EIR3", xcable4pm3, ycable4pm3, zcablePS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB3", 1,  "EIR7", xcable3pm3, ycable3pm3, zcableVeto, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB4", 1,  "EIR7", xcable4pm3, ycable4pm3, zcableVeto, 0, "ONLY");
  
  TVirtualMC::GetMC()->Gspos("ECB2", 1,  "EIR4", xcable21pm4,    0.,     zcablePS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB2", 2,  "EIR4", xcable22pm4,    0.,     zcablePS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB3", 1,  "EIR4", xcable3pm4, ycable3pm4, zcablePS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB2", 1,  "EIR8", xcable21pm4,    0.,     zcableVeto, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB2", 2,  "EIR8", xcable22pm4,    0.,     zcableVeto, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ECB3", 1,  "EIR8", xcable3pm4, ycable3pm4, zcableVeto, 0, "ONLY");
     


  //=============== NOW POSITIONING THE Al Boxes IN EPM'S================//
  
   
  TVirtualMC::GetMC()->Gspos("ESV1", 1,  "EPM1",  dboxFea1[0]  - esvdA1[0] - 8.0,  0., zelvdbPS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EVV1", 1,  "EPM1",  dboxFea1[0]  - esvdB1[0] - 8.0,  0., zelvdbVeto, 0, "ONLY");
  
  TVirtualMC::GetMC()->Gspos("ESV2", 1,  "EPM2", -dboxFea2[0]  + esvdA1[0] + 8.0, 2.3, zelvdbPS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EVV2", 1,  "EPM2", -dboxFea2[0]  + esvdB1[0] + 8.0, 2.3, zelvdbVeto, 0, "ONLY");
  
  TVirtualMC::GetMC()->Gspos("ESV3", 1,  "EPM3", -dboxFea3[0]  + esvdA1[0] + 8.0,  0., zelvdbPS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EVV3", 1,  "EPM3", -dboxFea3[0]  + esvdB1[0] + 8.0,  0., zelvdbVeto, 0, "ONLY");
  
  TVirtualMC::GetMC()->Gspos("ESV4", 1,  "EPM4",  dboxFea4[0]  - esvdA1[0] - 8.0, 2.3, zelvdbPS, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EVV4", 1,  "EPM4",  dboxFea4[0]  - esvdB1[0] - 8.0, 2.3, zelvdbVeto, 0, "ONLY");
  
  //==================================================================//
  //====================== LAST THING IS TO INSTALL ELMB ================//
  
  // ELMB,s are the G10 Volumes

  // First to create Air Volume to place ELMBs
  Float_t xelmb[3];
  xelmb[0] = 10.0;
  xelmb[1] = 4.0;
  xelmb[2] = 0.5;
  
  TVirtualMC::GetMC()->Gsvolu("ELMB", "BOX", idtmed[698], xelmb, 3);
  
  // There are more G10 Volumes
  // But in approximation, we reduced them to two
  // ELM1 & ELM2
  
  Float_t xelmb1[3];
  xelmb1[0] = 9.7;
  xelmb1[1] = 3.6;
  xelmb1[2] = 0.1;
  
  TVirtualMC::GetMC()->Gsvolu("ELM1", "BOX", idtmed[607], xelmb1, 3);
  
  Float_t xelmb2[3];
  xelmb2[0] = 6.0;
  xelmb2[1] = 3.0;
  xelmb2[2] = 0.1;
  
  TVirtualMC::GetMC()->Gsvolu("ELM2", "BOX", idtmed[607], xelmb2, 3);
  
  /******** NOW POSITIONING THE G10 VOLUMES ELM1 & ELM2 IN ELMB **********/
  
  TVirtualMC::GetMC()->Gspos("ELM1", 1,  "ELMB",  0., 0., -0.3, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ELM2", 1,  "ELMB",  0., 0.,  0.3, 0, "ONLY");
  
  // Position co-ordinates of ELMBs in EPM2 & EPM4 
  
  Float_t xelmbepm2, xelmbepm4, yelmbepm2, yelmbepm4, zelmbPS, zelmbVeto;
  
  xelmbepm2 = -gaspmd2[0] + 16.0 +23.2 + 2.5 + xelmb[0];
  xelmbepm4 =  gaspmd4[0] - 16.0 -23.2 - 2.5 - xelmb[0];
  
  yelmbepm2 = -gaspmd2[1] + 1.0 + xelmb[1];
  yelmbepm4 = -gaspmd4[1] + 1.0 + xelmb[1];
  
  zelmbPS   = zfea + fgkThSteel/2.+  xelmb[2];
  zelmbVeto = zfea - fgkThSteel/2.-  xelmb[2];
  
  /************ NOW PLACE ELMB'S IN EPM2 & EPM4 *********************/
  
  // There are total of 14 ELMB volumes
  // three on both sides of EPM2 (total of 6)
  // and four on both sides of EPM4 (total of 8)
  // The ELMBs are placed at the bottom of 
  // SS support, which is the extended part
  
  // Placement of ELMBs on EPM2
  for(Int_t kk=1;kk<=3;kk++){
    TVirtualMC::GetMC()->Gspos("ELMB", kk,  "EPM2",  xelmbepm2, yelmbepm2, zelmbPS, 0, "ONLY");
    xelmbepm2 = xelmbepm2 + xelmb[0] + 0.5 + xelmb[0];
  }
  
  xelmbepm2 = -gaspmd2[0] + 16.0 +23.2 + 2.5 + xelmb[0];
  
  for(Int_t kk=4;kk<=6;kk++){
    TVirtualMC::GetMC()->Gspos("ELMB", kk, "EPM2", xelmbepm2, yelmbepm2, zelmbVeto, 0, "ONLY");
    xelmbepm2 = xelmbepm2 + xelmb[0] + 0.5 + xelmb[0];
  }
  
  // Placement of ELMBs on EPM4
  for(Int_t kk=1;kk<=4;kk++){
    TVirtualMC::GetMC()->Gspos("ELMB", kk, "EPM4", xelmbepm4, yelmbepm4, zelmbPS, 0, "ONLY");
    xelmbepm4 = xelmbepm4 - xelmb[0] - 0.5 - xelmb[0];
  }
  
  xelmbepm4 =  gaspmd4[0] - 16.0 -23.2 - 2.5 - xelmb[0];
  for(Int_t kk=5;kk<=8;kk++){
    TVirtualMC::GetMC()->Gspos("ELMB", kk, "EPM4", xelmbepm4, yelmbepm4, zelmbVeto, 0, "ONLY");
    xelmbepm4 = xelmbepm4 - xelmb[0] - 0.5 - xelmb[0];
  }
  
  //========= Placement of ELMBs Completed ============================/
  
  // -------------  Now to Place Unit Modules in four quadrants 
  //                EPM1, EPM2, EPM3 & EPM4 ---------------------//

  // Position co-ordinates of Unit Modules
  
  Double_t xcord[24];
  Double_t ycord[24];
  
  xcord[0]  = xa1;
  xcord[1]  = xa2;
  xcord[2]  = xa3;
  xcord[3]  = xa1;
  xcord[4]  = xa2;
  xcord[5]  = xa3;
  xcord[6]  = -xa1;
  xcord[7]  = -xa2;
  xcord[8]  = -xa3;
  xcord[9]  = -xa1;
  xcord[10] = -xa2;
  xcord[11] = -xa3;
  xcord[12] = xb1;
  xcord[13] = xb2;
  xcord[14] = xb1;
  xcord[15] = xb2;
  xcord[16] = xb1;
  xcord[17] = xb2;
  xcord[18] = -xb1;
  xcord[19] = -xb2;
  xcord[20] = -xb1;
  xcord[21] = -xb2;
  xcord[22] = -xb1;
  xcord[23] = -xb2;

  ycord[0]  = ya1;
  ycord[1]  = ya1;
  ycord[2]  = ya1;
  ycord[3]  = ya2;
  ycord[4]  = ya2;
  ycord[5]  = ya2;
  ycord[6]  = -ya1;
  ycord[7]  = -ya1;
  ycord[8]  = -ya1;
  ycord[9]  = -ya2;
  ycord[10] = -ya2;
  ycord[11] = -ya2;
  ycord[12] = yb1;
  ycord[13] = yb1;
  ycord[14] = yb2;
  ycord[15] = yb2;
  ycord[16] = yb3+0.100007; //Because of overlapping the factor 0.100007 
  ycord[17] = yb3+0.100007; // is added
  ycord[18] = -yb1;
  ycord[19] = -yb1;
  ycord[20] = -yb2;
  ycord[21] = -yb2;
  ycord[22] = -yb3-0.100004; //Because of overlapping the factor 0.100007 
  ycord[23] = -yb3-0.100004; // is added
 

  // Placement of unit modules EUM1 & EUV1(long type)
  // and EUM2 & EUV2(short type)
  // in the four quadrants EPM1, EPM2, EPM3 & EPM4
  
  for(Int_t ii=0;ii<=5;ii++){
    if(fModStatus[ii]){
      TVirtualMC::GetMC()->Gspos("EUM1", ii, "EPM1", xcord[ii]+xLead1,ycord[ii]+yLead1, zpsa, 0, "ONLY");
    }  
  }
  
  for(Int_t ii=6;ii<=11;ii++){
    if(fModStatus[ii]) {
      TVirtualMC::GetMC()->Gspos("EUM1", ii, "EPM2", xcord[ii]+xLead2, ycord[ii]+yLead2, zpsa, jhrot12, "ONLY");
    }
  }
  
  for(Int_t ii=12;ii<=17;ii++){
    if(fModStatus[ii]) {
      TVirtualMC::GetMC()->Gspos("EUM2", ii, "EPM3", xcord[ii]+xLead3, ycord[ii]+yLead3, zpsa, 0, "ONLY");
    }
  }
  
  for(Int_t ii=18;ii<=23;ii++){
    if(fModStatus[ii]) {
      TVirtualMC::GetMC()->Gspos("EUM2", ii, "EPM4", xcord[ii]+xLead4, ycord[ii]+yLead4, zpsa, jhrot12, "ONLY");
    }
  }
  
  for(Int_t ii=24;ii<=29;ii++){
    if(fModStatus[ii]) {
      TVirtualMC::GetMC()->Gspos("EUV1", ii, "EPM1", xcord[ii-24]+xLead1, ycord[ii-24]+yLead1, zcva, 0, "ONLY");
    }
  }
  
  for(Int_t ii=30;ii<=35;ii++){
    if(fModStatus[ii]) {
      TVirtualMC::GetMC()->Gspos("EUV1", ii, "EPM2", xcord[ii-24]+xLead2, ycord[ii-24]+yLead2, zcva, jhrot12, "ONLY");
    }
  }
  
  for(Int_t ii=36;ii<=41;ii++){
    if(fModStatus[ii]) {
      TVirtualMC::GetMC()->Gspos("EUV2", ii, "EPM3", xcord[ii-24]+xLead3, ycord[ii-24]+yLead3, zcva, 0, "ONLY");
    }
  }
  
  for(Int_t ii=42;ii<=47;ii++){
    if(fModStatus[ii]) {
      TVirtualMC::GetMC()->Gspos("EUV2", ii, "EPM4", xcord[ii-24]+xLead4, ycord[ii-24]+yLead4, zcva, jhrot12, "ONLY");
    }
  }
  
  //-------------- Placement of Unit Modules Completed ---------------// 
  
  // ========== PLACE THE EPMD IN ALICE ======================//  
  
  // Now the Job to assemble the five mother volumes of PMD in ALICE
  
  // Z-distance of PMD from Interaction Point

  zp = fgkZdist;
  
  // X and Y-positions of the EPM1, EPM2, EPM3 & EPM4
  Float_t xfinal,yfinal; 
  Float_t xsm1,  xsm2,  xsm3,  xsm4;
  Float_t ysm1,  ysm2,  ysm3,  ysm4;
  
  xfinal = (fSMLengthax + serviceX/2. + serviceXext/2. + 0.05) + 0.48/2. +
    (fSMLengthbx + serviceX/2. + serviceXext/2.+ 0.05);

  //Extra width of the SS plate on Support Structure on X-side and 1mm thick SS for cooling encloser 
  //Extra width of the SS plate on Support Structure on X-side for B-Type
  
  yfinal = (fSMLengthay + serviceYa/2.)+ 0.20/2 + (fSMLengthby + serviceYb/2.);

  //serviceYa is the Extra width of the SS plate on Support Structur on Y-side for EPM1 & EPM3 
  //serviceYb is the Extra width of the SS plate on Support Structur on Y-side for EPM2 & EPM4
  
  
  xsm1 =  xfinal  - (fSMLengthax + serviceX/2. + serviceXext/2. + 0.05);
  ysm1 =  yfinal  - (fSMLengthay + serviceYa/2.) - 2.3;
  
  xsm2 =  -xfinal  + (fSMLengthax + serviceX/2. + serviceXext/2. + 0.05);
  ysm2 =  -yfinal  + (fSMLengthay + serviceYb/2.) - 2.3;
  
  xsm3 =  -xfinal + (fSMLengthbx + serviceX/2. + serviceXext/2. + 0.05);
  ysm3 =   yfinal - (fSMLengthby + serviceYa/2.) - 2.3;
  
  xsm4 =   xfinal - (fSMLengthbx + serviceX/2. + serviceXext/2. + 0.05);
  ysm4 =  -yfinal + (fSMLengthby + serviceYb/2.) - 2.3;
  
  //Position Full PMD in ALICE   
  //
  //       EPM1                EPM3
  //
  //       EPM4                EPM2
  //  (rotated EPM3)      (rotated EPM1)
  //
  //                EFGD
  //        (Girders and its Carriage)
  
  TVirtualMC::GetMC()->Gspos("EPM1", 1, "ALIC",  xsm1,ysm1,zp, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPM2", 1, "ALIC",  xsm2,ysm2,zp, 0, "ONLY");
 
  TVirtualMC::GetMC()->Gspos("EPM3", 1, "ALIC",  xsm3,ysm3,zp, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("EPM4", 1, "ALIC",  xsm4,ysm4,zp, 0, "ONLY");
  
  TVirtualMC::GetMC()->Gspos("EFGD", 1, "ALIC", 0., yfinal + fulgrdr[1], zp, 0, "ONLY");  
}

//_____________________________________________________________________________

void AliPMDv1::CreateMaterials()
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
  
  
  // --- CH2 : PLASTIC  --- 
  
  Float_t aCH2[2] = { 12.,1.};
  Float_t zCH2[2] = { 6.,1.};
  Float_t wCH2[2] = { 1.,2.};
  Float_t dCH2    = 0.95;
  AliMixture(31, "CH2  $", aCH2, zCH2, dCH2, -2, wCH2);
  
  // --- CABLES : 80% Plastic and 20% Copper  --- 
  
  Float_t aCABLE[3] = { 12.,1.,63.5 };
  Float_t zCABLE[3] = { 6.,1.,29. };
  Float_t wCABLE[3] = { 0.6857, 0.1143, 0.2};
  Float_t dCABLE    = dCH2*0.8 + 8.96*0.2;
  AliMixture(32, "CABLE  $", aCABLE, zCABLE, dCABLE, 3, wCABLE);

  
  
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
  AliMedium(32, "CABLE   $", 32, 0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
  AliMedium(98, "Vacuum  $", 98, 0, 0, isxfld, sxmgmx, 1., .1, .10, 10);
  AliMedium(99, "Air gaps$", 99, 0, 0, isxfld, sxmgmx, 1., .1, .10, .1);
  
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
  // --- Generate explicitly delta rays in the iron, aluminium and lead --- 
  // Gstpar is removed from this place and 
  // the energy cut offs in the medium moved to galice.cuts
  
  //TVirtualMC::GetMC()->Gstpar(idtmed[605], "LOSS", 3.);
  //TVirtualMC::GetMC()->Gstpar(idtmed[605], "DRAY", 1.);
  
  // Visualization of volumes
  gGeoManager->SetVolumeAttribute("ECAR", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECCU", "SEEN", 1);
  gGeoManager->SetVolumeAttribute("ECCU", "COLO", 4);
  gGeoManager->SetVolumeAttribute("EST1", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EST2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EHC1", "SEEN", 0);  
  gGeoManager->SetVolumeAttribute("EHC2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EDGA", "SEEN", 1);
  gGeoManager->SetVolumeAttribute("EDGB", "SEEN", 1);
  gGeoManager->SetVolumeAttribute("EEGA", "SEEN", 1);
  gGeoManager->SetVolumeAttribute("EEGB", "SEEN", 1);
  gGeoManager->SetVolumeAttribute("EUM1", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EUV1", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EUM2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EUV2", "SEEN", 0);

 
  gGeoManager->SetVolumeAttribute("EFEE", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EFEE", "COLO", 4);
  gGeoManager->SetVolumeAttribute("EFBA", "SEEN", 1);
  gGeoManager->SetVolumeAttribute("EFBA", "COLO", 4);
  gGeoManager->SetVolumeAttribute("EFBB", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EFBB", "COLO", 4);

  gGeoManager->SetVolumeAttribute("ELDA", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ELDB", "SEEN", 0);

  gGeoManager->SetVolumeAttribute("EFE1", "SEEN", 0); 
  gGeoManager->SetVolumeAttribute("EFE2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EFE3", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EFE4", "SEEN", 0);

  gGeoManager->SetVolumeAttribute("ESC1", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECC1", "COLO", 2);
  gGeoManager->SetVolumeAttribute("ESC2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECC2", "COLO", 2);
  gGeoManager->SetVolumeAttribute("ESC3", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECC3", "COLO", 2);
  gGeoManager->SetVolumeAttribute("ESC4", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECC4", "COLO", 2);

  gGeoManager->SetVolumeAttribute("ECC1", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECC2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECC3", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECC4", "SEEN", 0);

  gGeoManager->SetVolumeAttribute("EPM1", "SEEN", 1);
  gGeoManager->SetVolumeAttribute("EPM2", "SEEN", 1);
  gGeoManager->SetVolumeAttribute("EPM3", "SEEN", 1);
  gGeoManager->SetVolumeAttribute("EPM4", "SEEN", 1);

  gGeoManager->SetVolumeAttribute("ECB1", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECB2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECB3", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ECB4", "SEEN", 0);

  gGeoManager->SetVolumeAttribute("ELMB", "SEEN", 0);
  
  gGeoManager->SetVolumeAttribute("ESV1", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ESV2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ESV3", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("ESV4", "SEEN", 0);

  gGeoManager->SetVolumeAttribute("EVV1", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EVV2", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EVV3", "SEEN", 0);
  gGeoManager->SetVolumeAttribute("EVV4", "SEEN", 0);

  gGeoManager->SetVolumeAttribute("EFGD", "SEEN", 0);
}

//_____________________________________________________________________________

void AliPMDv1::StepManager()
{
  //
  // Called at each step in the PMD
  //
  
  Int_t   copy;
  Float_t hits[5], destep;
  Float_t center[3] = {0,0,0};
  Int_t   vol[6];
  //const char *namep;
  //    printf("Current vol  is ********  %s \n",namep);
  if(fMC->CurrentMedium() == fMedSens && (destep = fMC->Edep())) {
    
    fMC->CurrentVolID(copy);
    //namep=fMC->CurrentVolName();
    //  printf("Current vol  is %s \n",namep);
    vol[0]=copy;
    
    fMC->CurrentVolOffID(1,copy);
    //namep=fMC->CurrentVolOffName(1);
    // printf("Current vol 11 is %s \n",namep);
    vol[1]=copy;
    
    fMC->CurrentVolOffID(2,copy);
    //namep=fMC->CurrentVolOffName(2);
    // printf("Current vol 22 is %s \n",namep);
    vol[2]=copy;
    
    fMC->CurrentVolOffID(3,copy);
    //namep=fMC->CurrentVolOffName(3);
    // printf("Current vol 33 is %s \n",namep);
    vol[3]=copy;
    
    fMC->CurrentVolOffID(4,copy);
    //namep=fMC->CurrentVolOffName(4);
    // printf("Current vol 44 is %s \n",namep);
    vol[4]=copy;
    
    fMC->CurrentVolOffID(5,copy);
    //namep=fMC->CurrentVolOffName(5);
    //printf("Current vol 55 is %s \n",namep);
    vol[5]=copy;

    
    // printf("volume number %4d %4d %4d %4d %4d %4d %10.3f \n",vol[0],vol[1],vol[2],vol[3],vol[4],vol[5],destep*1000000);// edep in MeV
    
    
    fMC->Gdtom(center,hits,1);
    hits[3] = destep*1e9; //Number in eV

    // this is for pile-up events
    hits[4] = fMC->TrackTime();

    AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);

    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kPMD);

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
  
  fSMLengthbx = 42.6136;
  //The total length in X is due to the following components
  // Factor 2 is because of 2 module length in X for this type
  // fgkNcolUM2*fgkCellRadius (96 x 0.25): Total span of each module in X
  // fgkCellRadius/2. : There is offset of 1/2 cell
  // 0.05+0.05 : Insulation gaps etc
  // fgkSSBoundary (0.3) : Boundary frame
  //double XB = 2.0*((fgkCellRadius/fgkSqroot3by2*fgkNcolUM2)-(fgkCellRadius*fgkSqroot3*(fgkNcolUM2-1)/6.)+(2.0*fgkGap)+(2.0*fgkGap)+fgkSSBoundary) + 0.1; 
  

  
  fSMLengthay = 49.35;
  //The total length in Y is due to the following components
  // Factor 2 is because of 2 module length in Y for this type
  // fgkCellRadius/fgkSqroot3by2)*fgkNrowUM1 (0.25/sqrt3/2 * 96): Total span of each module in Y
  //  of strips
  // 0.05+0.05 : Insulation gaps etc
  // fgkSSBoundary (0.3) : Boundary frame
  // 0.6cm is the channel width plus tolerance
  // double  YA = 2.0*(fgkNrowUM1*fgkCellRadius+fgkCellRadius/2.+(2.0*fgkGap)+(2.0*fgkGap)+fgkSSBoundary) +  0.6/2.;
  
  fSMLengthby =  37.925;
  //The total length in Y is due to the following components
  // Factor 3 is because of 3 module length in Y for this type
  // fgkCellRadius/fgkSqroot3by2)*fgkNrowUM2 (0.25/sqrt3/2 * 48): Total span of each module in Y
  //  of strips
  // 0.05+0.05 : Insulation gaps etc
  // fgkSSBoundary (0.3) : Boundary frame
  // 10mm is the channel width plus tolerance
  //double YB = 3.0*((fgkNrowUM2*fgkCellRadius + fgkCellRadius/2.)+(2.0*fgkGap)+(2.0*fgkGap)+fgkSSBoundary) + 1.0/2.;
  
  
  //Thickness of a pre/veto plane 
  fDthick     = fgkThSS/2. + 1.2;     // 1.2 added as FEE Board are now assembled with pre/veto
  
  //Thickness of the PMD ; 2.4 added for FEE boards 
  fSMthickpmd    = 2.0*(fgkThSS/2.) +fgkThSteel/2.+fgkThLead/2.0 + 2.4/2.;
  
  fSMthick = 17.; //17 cm is the full profile of PMD
  
  
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
    //for(Int_t cnt=1; cnt<=4; cnt++){
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
void AliPMDv1::SetCpvOff()
{
  // Set the entire CPV plane off

  for (Int_t imodule = 24; imodule < 48; imodule++)
    fModStatus[imodule] = 0;
}
// ------------------------------------------------------------------
void AliPMDv1::SetPreOff()
{
  // Set the entire Preshower plane off

  for (Int_t imodule = 0; imodule < 24; imodule++)
    fModStatus[imodule] = 0;

}
// ------------------------------------------------------------------
void AliPMDv1::SetModuleOff(Int_t imodule)
{
  // Set the individual module off

  fModStatus[imodule] = 0;

}
