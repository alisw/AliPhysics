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

//////////////////////////////////////////////////////////////////////
//                                                                  //
//  (V-zero) detector  version 0  as designed by the Lyon group     //
//   All comments should be sent to Brigitte CHEYNIS :              //
//                                  b.cheynis@ipnl.in2p3.fr         // 
//   Geometrie      du 25/02/2002                                   //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <string.h>
#include <Riostream.h>

#include <TBRIK.h>
#include <TBox.h>
#include <TCONE.h>
#include <TClonesArray.h>
#include <TGeometry.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TNode.h>
#include <TObjectTable.h>
#include <TPCON.h>
#include <TPGON.h>
#include <TSPHE.h>
#include <TShape.h>
#include <TTRAP.h>
#include <TTRD2.h>
#include <TVirtualMC.h>

#include "ABSOConst.h"
#include "ABSOSHILConst.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliVZEROdigit.h"
#include "AliVZEROhit.h"
#include "AliVZEROv0.h"
#include "AliMC.h"

ClassImp(AliVZEROv0)

//--------------------------------------------------------------------
AliVZEROv0:: AliVZEROv0():AliVZERO()
{


}
//--------------------------------------------------------------------
AliVZEROv0::AliVZEROv0(const char *name, const char *title):
 AliVZERO(name,title)
{

// Standard constructor for V-zeroR Detector (right part)  version 0


  Int_t i;

  printf("\n");
  for(i=0;i<30;i++) printf("*");
  printf(" Create VZERO object ");
  for(i=0;i<30;i++) printf("*");
  printf("\n");
  
}

//-------------------------------------------------------------------------
void AliVZEROv0::CreateGeometry()
{

// Creates the Geant geometry of the V-zero Detector  version 0

  
  Int_t i;
  
  printf("\n");
  for(i=0;i<30;i++) printf("*");
  printf(" Create VZERO Geometry ");
  for(i=0;i<30;i++) printf("*");
  printf("\n");
    
  
  Int_t    *idtmed = fIdtmed->GetArray()-2999;

  Int_t    n_detec_R = 1;
  Int_t    n_detec_L = 1;
 
  Int_t    n_cells_R = 1;
  Int_t    n_cells_L = 1;
  
  Int_t    idrotm[999];
 
  Float_t  height1, height2, height3, height4, height5; 
  Float_t  height;
  Float_t  theta;  
  Float_t  half_thick_alu;
  Float_t  half_thick_qua1,half_thick_qua2,half_thick_qua3;
  Float_t  half_thick_qua4,half_thick_qua5;
  Float_t  zdet;
  Float_t  r0, r5;
  Float_t  pi = TMath::Pi();
  Float_t  thick_alu;
  
  height1           =     2.0;           // height of cell 1, in cm
  height2           =     3.2;           // height of cell 2, in cm
  height3           =     4.9;           // height of cell 3, in cm
  height4           =     7.5;           // height of cell 4, in cm
  height5           =    12.0;           // height of cell 5, in cm
  
  theta             = pi/6.0/2.0;       // half angular opening = 15 degrees
  half_thick_alu    = 0.0025;            // half thickness of aluminum foil, in cm
  thick_alu         = 2.0 * half_thick_alu; 
  fThickness1	    = 2.5;
  half_thick_qua1   = fThickness1/2.0;   // half thickness of WRAPPED quartz cell (inner ring)
  half_thick_qua2   = half_thick_qua1  - 0.25;
  half_thick_qua3   = half_thick_qua2  - 0.25;
  half_thick_qua4   = half_thick_qua3  - 0.25;
  half_thick_qua5   = half_thick_qua4  - 0.25;
  
  zdet              =    86.9 +fThickness/2.0;  // distance to vertex (along Z axis)
  r0                =    4.0;            // closest distance to center of the beam pipe
  height            =    height1 + height2 + height3 + height4 + height5;
  r5                =    r0 + height;

//............................................................................

// Here I add the flange which is sitting on beam line 
// right in front of V0R detector, and which I found on CERN drawing 
// entitled : ALICE BEAM VACCUM CHAMBER - RB26 version III :  
   
//     Float_t   pflange[3];
//     
//     pflange[0] = 3.0;
//     pflange[1] = 5.675;
//     pflange[2] = 0.9;        
// 
//     gMC->Gsvolu("QFA0","TUBE", idtmed[3003], pflange, 3);
//     gMC->Gspos("QFA0", 1 ,"ALIC", 0.0, 0.0, 85.0+0.9, 0, "ONLY");
     
//............................................................................


// Creation of mother volume V0LE - left part - :
// Face entree a -350.0 cm ...

   Float_t   partube[3];
   
   partube[0] =  4.3;
   partube[1] = 45.0;
   partube[2] = fThickness1/2.0;   
    
   gMC->Gsvolu("V0LE","TUBE",idtmed[3002],partube,3);
  
   
// Creation of five rings - left part - :
// Face entree a -350.0 cm ... 

// Mother volume V0L0 in which will be set 5 quartz cells 


  Float_t   par[11];
  
  Float_t   dist0_left;   
  Float_t   r0_left      =   4.3;   
  Float_t   height1_left =   2.6; 
  Float_t   height2_left =   4.1;
  Float_t   height3_left =   6.4;
  Float_t   height4_left =  10.2;
  Float_t   height5_left =  16.9;
  Float_t   height_left  = height1_left + height2_left + height3_left 
                                        + height4_left + height5_left;
  Float_t   r5_left      = r0_left  + height_left; 
  
   
  dist0_left  =  r0_left + height_left / 2.0;
  thick_alu   =  2.0*half_thick_alu;
  
  par[0]      =  half_thick_qua1;
  par[1]      =  0.0;
  par[2]      =  0.0;
  par[3]      =  height_left / 2.0 ;
  par[4]      =  TMath::Tan(theta) * r0_left;
  par[5]      =  TMath::Tan(theta) * r5_left;
  par[6]      =  0.0;
  par[7]      =  height_left / 2.0 ;
  par[8]      =  TMath::Tan(theta) * r0_left;
  par[9]      =  TMath::Tan(theta) * r5_left;
  par[10]     =  0.0;
  

  gMC->Gsvolu("V0L0","TRAP",idtmed[3010],par,11);  // air volume
  
  Float_t   dist1_left;
  Float_t   r1_left;
  Float_t   offset_left;
     
  dist1_left     =  (- height_left + height1_left) /2.0; 
  r1_left        =  r0_left + height1_left;
  offset_left    = - fThickness1/2.0 + 0.1; 
   
  par[0]    =  half_thick_qua1 - thick_alu;
  par[3]    =  height1_left / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r0_left - thick_alu;
  par[5]    =  TMath::Tan(theta) * r1_left- thick_alu;
  par[7]    =  height1_left / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r0_left - thick_alu;
  par[9]    =  TMath::Tan(theta) * r1_left - thick_alu;


  gMC->Gsvolu("V0L1","TRAP",idtmed[3002],par,11);  // quartz volume
  gMC->Gspos("V0L1",1,"V0L0", 0.0, dist1_left , 0.0, 0,"ONLY"); 
  
  Float_t   dist2_left;
  Float_t   r2_left; 
    
  dist2_left     =    (- height_left + height2_left) /2.0 + height1_left;
  r2_left        =       r1_left + height2_left; 
  
  par[0]    =  half_thick_qua1 - thick_alu;   
  par[3]    =  height2_left / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r1_left - thick_alu;
  par[5]    =  TMath::Tan(theta) * r2_left - thick_alu;
  par[7]    =  height2_left / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r1_left - thick_alu;
  par[9]    =  TMath::Tan(theta) * r2_left - thick_alu;

  gMC->Gsvolu("V0L2","TRAP",idtmed[3002],par,11);  // quartz volume
  gMC->Gspos("V0L2",1,"V0L0", 0.0, dist2_left , 0.0, 0,"ONLY"); 
  
  
  Float_t   dist3_left;
  Float_t   r3_left;
     
  dist3_left     =    (- height_left + height3_left) /2.0 + height1_left + height2_left;
  r3_left        =       r2_left + height3_left; 
   
  par[0]    =  half_thick_qua1 - thick_alu;   
  par[3]    =  height3_left / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r2_left - thick_alu;
  par[5]    =  TMath::Tan(theta) * r3_left - thick_alu;
  par[7]    =  height3_left / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r2_left - thick_alu;
  par[9]    =  TMath::Tan(theta) * r3_left - thick_alu;

  gMC->Gsvolu("V0L3","TRAP",idtmed[3002],par,11);  // quartz volume
  gMC->Gspos("V0L3",1,"V0L0", 0.0, dist3_left , 0.0, 0,"ONLY");
    
  Float_t   dist4_left;
  Float_t   r4_left;
       
  dist4_left     =    (- height_left + height4_left) /2.0 + height1_left 
                                     + height2_left + height3_left;
  r4_left        =       r3_left + height4_left; 
   
  par[0]    =  half_thick_qua1 - thick_alu;   
  par[3]    =  height4_left / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r3_left - thick_alu;
  par[5]    =  TMath::Tan(theta) * r4_left - thick_alu;
  par[7]    =  height4_left / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r3_left - thick_alu;
  par[9]    =  TMath::Tan(theta) * r4_left - thick_alu;

  gMC->Gsvolu("V0L4","TRAP",idtmed[3002],par,11);  // quartz volume
  gMC->Gspos("V0L4",1,"V0L0", 0.0, dist4_left , 0.0, 0,"ONLY");


  Float_t   dist5_left;

         
  dist5_left     =    (- height_left + height5_left) /2.0 + height1_left 
                                     + height2_left + height3_left + height4_left;

   
  par[0]    =  half_thick_qua1 - thick_alu;   
  par[3]    =  height5_left / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r4_left - thick_alu;
  par[5]    =  TMath::Tan(theta) * r5_left - thick_alu;
  par[7]    =  height5_left / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r4_left - thick_alu;
  par[9]    =  TMath::Tan(theta) * r5_left - thick_alu;

  gMC->Gsvolu("V0L5","TRAP",idtmed[3002],par,11);  // quartz volume
  gMC->Gspos("V0L5",1,"V0L0", 0.0, dist5_left , 0.0, 0,"ONLY");
  
  
//............................................................................

// Creation of mother volume V0RI - right part - :

  
  partube[0] = r0 - 0.2;
  partube[1] = (r5 + 1.0) / TMath::Cos(theta);
  partube[2] = fThickness/2.0; 
  
    
  gMC->Gsvolu("V0RI","TUBE",idtmed[3010],partube,3);
  
// Creation of  carbon lids (1 mm thick) to keep V0RI box shut...

  Float_t   parbox[10];
  
  parbox[0] =    0.;
  parbox[1] =  360.;
  parbox[2] =    12;
  parbox[3] =     2;
  parbox[4] =  -0.1/2.0;
  parbox[5] =  r0;
  parbox[6] =  r5;     
  parbox[7] =  +0.1/2.0;
  parbox[8] =  r0;
  parbox[9] =  r5;  
  
  
  gMC->Gsvolu("V0CA","PGON",idtmed[3001],parbox,10); 
  gMC->Gspos("V0CA",1,"V0RI",0.0,0.0, fThickness/2.0-parbox[7],0,"ONLY");
  gMC->Gspos("V0CA",2,"V0RI",0.0,0.0,-fThickness/2.0+parbox[7],0,"ONLY");
  
// Creation of aluminum rings to maintain the V0RI pieces ...

  parbox[4] =  -fThickness/2.0;
  parbox[5] =  r0 -0.2;
  parbox[6] =  r0;     
  parbox[7] =  +fThickness/2.0;
  parbox[8] =  r0 -0.2;
  parbox[9] =  r0; 
  
  gMC->Gsvolu("V0IR","PGON",idtmed[3003],parbox,10);    
  gMC->Gspos("V0IR",1,"V0RI",0.0,0.0,0.0,0,"ONLY");
  
  parbox[4] =  -fThickness/2.0;
  parbox[5] =  r5;
  parbox[6] =  r5 + 1.0;     
  parbox[7] =  +fThickness/2.0;
  parbox[8] =  r5;
  parbox[9] =  r5 + 1.0; 
 
  gMC->Gsvolu("V0ER","PGON",idtmed[3003],parbox,10);    
  gMC->Gspos("V0ER",1,"V0RI",0.0,0.0,0.0,0,"ONLY");
  
// Mother volume V0R0 in which will be set 5  quartz cells 
// each one  WRAPPED in reflecting aluminum : 
  
  Float_t   dist0;   
    
  dist0     =  r0 + height / 2.0;
  thick_alu =  2.0*half_thick_alu;
  
  par[0]    =  half_thick_qua1;
  par[1]    =  0.0;
  par[2]    =  0.0;
  par[3]    =  height / 2.0 ;
  par[4]    =  TMath::Tan(theta) * r0;
  par[5]    =  TMath::Tan(theta) * r5;
  par[6]    =  0.0;
  par[7]    =  height / 2.0 ;
  par[8]    =  TMath::Tan(theta) * r0;
  par[9]    =  TMath::Tan(theta) * r5;
  par[10]   =  0.0;


  gMC->Gsvolu("V0R0","TRAP",idtmed[3010],par,11);  // air volume

// Elementary cell of ring 1 :
 
  Float_t   dist1;
  Float_t   r1;
  Float_t   offset;
     
  dist1     =  (- height + height1) /2.0; 
  r1        =  r0 + height1;
  offset    = - fThickness/2.0 + 0.1; 
   
  par[0]    =  half_thick_qua1 - thick_alu;
  par[3]    =  height1 / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r0 - thick_alu;
  par[5]    =  TMath::Tan(theta) * r1- thick_alu;
  par[7]    =  height1 / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r0 - thick_alu;
  par[9]    =  TMath::Tan(theta) * r1 - thick_alu;


  gMC->Gsvolu("V0R1","TRAP",idtmed[3002],par,11);  // quartz volume
  gMC->Gspos("V0R1",1,"V0R0", 0.0, dist1 , 0.0, 0,"ONLY");
  
  par[0]    =  half_thick_alu;

  gMC->Gsvolu("V0A1","TRAP",idtmed[3004],par,11);  // aluminum trap-shaped foil
  gMC->Gspos("V0A1",1,"V0R1",0.0,0.0, - half_thick_qua1 + half_thick_alu,0,"ONLY");  
  gMC->Gspos("V0A1",2,"V0R1",0.0,0.0, + half_thick_qua1 - half_thick_alu,0,"ONLY");  

  parbox[0] = half_thick_alu;
  parbox[1] = height1 / TMath::Cos(theta)/ 2.0;
  parbox[2] = half_thick_qua1;
  
  gMC->Gsvolu("V0A2","BOX",idtmed[3004],parbox,3);   // aluminum rectangular foil
  Float_t  theta_deg = 180./6./2.0; 
  Float_t h1;
  h1 = TMath::Tan(theta) * (r0 + height1/2.0);  
  AliMatrix(idrotm[911],90.0,+theta_deg,90.0,90.+theta_deg,0.0,0.);  
  gMC->Gspos("V0A2",1,"V0R1",-h1 + half_thick_alu,0.0,0.0,idrotm[911],"ONLY"); 
  AliMatrix(idrotm[912],90.0,-theta_deg,90.0,90.-theta_deg,0.0,0.);  
  gMC->Gspos("V0A2",2,"V0R1",+h1 - half_thick_alu,0.0,0.0,idrotm[912],"ONLY"); 

  parbox[0] = TMath::Tan(theta) * r0;
  parbox[1] = half_thick_alu;
  parbox[2] = half_thick_qua1;
  gMC->Gsvolu("V0A3","BOX",idtmed[3004],parbox,3);
  gMC->Gspos("V0A3",1,"V0R1",0.0, - (height1/2.0) + half_thick_alu ,0.0,0,"ONLY");
    

  parbox[0] = TMath::Tan(theta) * (r0 + height1);
  parbox[1] = half_thick_alu;
  parbox[2] = half_thick_qua1;
  gMC->Gsvolu("V0A4","BOX",idtmed[3004],parbox,3);
  gMC->Gspos("V0A4",1,"V0R1",0.0, (height1/2.0) - half_thick_alu,0.0,0,"ONLY");

  
//  Elementary cell of ring 2 : 
 
  Float_t   dist2;
  Float_t   r2; 
    
  dist2     =    (- height + height2) /2.0 + height1;
  r2        =       r1 + height2; 
  
  par[0]    =  half_thick_qua2 - thick_alu;   
  par[3]    =  height2 / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r1 - thick_alu;
  par[5]    =  TMath::Tan(theta) * r2 - thick_alu;
  par[7]    =  height2 / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r1 - thick_alu;
  par[9]    =  TMath::Tan(theta) * r2 - thick_alu;

  gMC->Gsvolu("V0R2","TRAP",idtmed[3002],par,11);  // quartz volume
  gMC->Gspos("V0R2",1,"V0R0", 0.0, dist2 , - half_thick_qua1 + half_thick_qua2, 0,"ONLY");

  par[0]    =  half_thick_alu;

  gMC->Gsvolu("V0B1","TRAP",idtmed[3004],par,11);  // aluminum trap-shaped foil
  gMC->Gspos("V0B1",1,"V0R2",0.0,0.0, - half_thick_qua2 + half_thick_alu,0,"ONLY");  
  gMC->Gspos("V0B1",2,"V0R2",0.0,0.0, + half_thick_qua2 - half_thick_alu,0,"ONLY");  

  parbox[0] = half_thick_alu;
  parbox[1] = height2 / TMath::Cos(theta)/ 2.0;
  parbox[2] = half_thick_qua2;
  
  gMC->Gsvolu("V0B2","BOX",idtmed[3004],parbox,3);   // aluminum rectangular foil
  Float_t h2;
  h2 = TMath::Tan(theta) * (r0  + height1 + height2/2.0);  
  gMC->Gspos("V0B2",1,"V0R2",-h2 + half_thick_alu,0.0,0.0,idrotm[911],"ONLY"); 
  gMC->Gspos("V0B2",2,"V0R2",+h2 - half_thick_alu,0.0,0.0,idrotm[912],"ONLY"); 

  parbox[0] = TMath::Tan(theta) * (r0 + height1);
  parbox[1] = half_thick_alu;
  parbox[2] = half_thick_qua2;
  gMC->Gsvolu("V0B3","BOX",idtmed[3004],parbox,3);
  gMC->Gspos("V0B3",1,"V0R2",0.0, - (height2/2.0) + half_thick_alu ,0.0,0,"ONLY");
    

  parbox[0] = TMath::Tan(theta) * (r0 + height1 +  height2);
  parbox[1] = half_thick_alu;
  parbox[2] = half_thick_qua2;
  gMC->Gsvolu("V0B4","BOX",idtmed[3004],parbox,3);
  gMC->Gspos("V0B4",1,"V0R2",0.0, (height2/2.0) - half_thick_alu,0.0,0,"ONLY");
  

// Elementary cell  ring 3 :   
  
  Float_t   dist3;
  Float_t   r3; 
    
  dist3     =    (- height + height3) /2.0 + height1 + height2;
  r3        =       r2 + height3; 
   
  par[0]    =  half_thick_qua3 - thick_alu;   
  par[3]    =  height3 / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r2 - thick_alu;
  par[5]    =  TMath::Tan(theta) * r3 - thick_alu;
  par[7]    =  height3 / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r2 - thick_alu;
  par[9]    =  TMath::Tan(theta) * r3 - thick_alu;

  gMC->Gsvolu("V0R3","TRAP",idtmed[3002],par,11);  // quartz volume
  gMC->Gspos("V0R3",1,"V0R0", 0.0, dist3 , - half_thick_qua1 +  half_thick_qua3, 0,"ONLY");


  par[0]    =  half_thick_alu;

  gMC->Gsvolu("V0C1","TRAP",idtmed[3004],par,11);  // aluminum trap-shaped foil
  gMC->Gspos("V0C1",1,"V0R3",0.0,0.0, - half_thick_qua3 + half_thick_alu,0,"ONLY");  
  gMC->Gspos("V0C1",2,"V0R3",0.0,0.0, + half_thick_qua3 - half_thick_alu,0,"ONLY");  

  parbox[0] = half_thick_alu;
  parbox[1] = height3 / TMath::Cos(theta)/ 2.0;
  parbox[2] = half_thick_qua3;
  
  gMC->Gsvolu("V0C2","BOX",idtmed[3004],parbox,3);   // aluminum rectangular foil
  Float_t h3;
  h3 = TMath::Tan(theta) * (r0  + height1 + height2 + height3/2.0);  
  gMC->Gspos("V0C2",1,"V0R3",-h3 + half_thick_alu,0.0,0.0,idrotm[911],"ONLY"); 
  gMC->Gspos("V0C2",2,"V0R3",+h3 - half_thick_alu,0.0,0.0,idrotm[912],"ONLY"); 

  parbox[0] = TMath::Tan(theta) * (r0 + height1 + height2);
  parbox[1] = half_thick_alu;
  parbox[2] = half_thick_qua3;
  gMC->Gsvolu("V0C3","BOX",idtmed[3004],parbox,3);
  gMC->Gspos("V0C3",1,"V0R3",0.0, - (height3/2.0) + half_thick_alu ,0.0,0,"ONLY");
    

  parbox[0] = TMath::Tan(theta) * (r0 + height1 + height2 + height3);
  parbox[1] = half_thick_alu;
  parbox[2] = half_thick_qua3;
  gMC->Gsvolu("V0C4","BOX",idtmed[3004],parbox,3);
  gMC->Gspos("V0C4",1,"V0R3",0.0, (height3/2.0) - half_thick_alu,0.0,0,"ONLY");


// Elementary cell  ring 4 :  

  Float_t   dist4;
  Float_t   r4; 
    
  dist4     =    (- height + height4) /2.0 + height1 + height2 + height3;
  r4        =       r3 + height4; 
   
  par[0]    =  half_thick_qua4 - thick_alu;   
  par[3]    =  height4 / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r3 - thick_alu;
  par[5]    =  TMath::Tan(theta) * r4 - thick_alu;
  par[7]    =  height4 / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r3 - thick_alu;
  par[9]    =  TMath::Tan(theta) * r4 - thick_alu;

  gMC->Gsvolu("V0R4","TRAP",idtmed[3002],par,11);  // quartz volume
  gMC->Gspos("V0R4",1,"V0R0", 0.0, dist4 , - half_thick_qua1 + half_thick_qua4, 0,"ONLY"); 


  par[0]    =  half_thick_alu;

  gMC->Gsvolu("V0D1","TRAP",idtmed[3004],par,11);  // aluminum trap-shaped foil
  gMC->Gspos("V0D1",1,"V0R4",0.0,0.0, - half_thick_qua4 + half_thick_alu,0,"ONLY");  
  gMC->Gspos("V0D1",2,"V0R4",0.0,0.0, + half_thick_qua4 - half_thick_alu,0,"ONLY");  

  parbox[0] = half_thick_alu;
  parbox[1] = height4 / TMath::Cos(theta)/ 2.0;
  parbox[2] = half_thick_qua4;
  
  gMC->Gsvolu("V0D2","BOX",idtmed[3004],parbox,3);   // aluminum rectangular foil
  Float_t h4;
  h4 = TMath::Tan(theta) * (r0  + height1 + height2 + height3 + height4/2.0);  
  gMC->Gspos("V0D2",1,"V0R4",-h4 + half_thick_alu,0.0,0.0,idrotm[911],"ONLY"); 
  gMC->Gspos("V0D2",2,"V0R4",+h4 - half_thick_alu,0.0,0.0,idrotm[912],"ONLY"); 

  parbox[0] = TMath::Tan(theta) * (r0 + height1 + height2 + height3);
  parbox[1] = half_thick_alu;
  parbox[2] = half_thick_qua4;
  gMC->Gsvolu("V0D3","BOX",idtmed[3004],parbox,3);
  gMC->Gspos("V0D3",1,"V0R4",0.0, - (height4/2.0) + half_thick_alu ,0.0,0,"ONLY");
    

  parbox[0] = TMath::Tan(theta) * (r0 + height1 + height2 + height3 + height4);
  parbox[1] = half_thick_alu;
  parbox[2] = half_thick_qua4;
  gMC->Gsvolu("V0D4","BOX",idtmed[3004],parbox,3);
  gMC->Gspos("V0D4",1,"V0R4",0.0, (height4/2.0) - half_thick_alu,0.0,0,"ONLY");


// Elementary cell  ring 5 :  

  Float_t   dist5;
      
  dist5     =    (- height + height5) /2.0 + height1 + height2 + height3 + height4;
  
  par[0]    =  half_thick_qua5 - thick_alu;      
  par[3]    =  height5 / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r4 - thick_alu;
  par[5]    =  TMath::Tan(theta) * r5 - thick_alu;
  par[7]    =  height5 / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r4 - thick_alu;
  par[9]    =  TMath::Tan(theta) * r5 - thick_alu;

  gMC->Gsvolu("V0R5","TRAP",idtmed[3002],par,11);  // quartz volume
  gMC->Gspos("V0R5",1,"V0R0", 0.0, dist5 , - half_thick_qua1 +  half_thick_qua5, 0,"ONLY"); 


  par[0]    =  half_thick_alu;

  gMC->Gsvolu("V0E1","TRAP",idtmed[3004],par,11);  // aluminum trap-shaped foil
  gMC->Gspos("V0E1",1,"V0R5",0.0,0.0, - half_thick_qua5 + half_thick_alu,0,"ONLY");  
  gMC->Gspos("V0E1",2,"V0R5",0.0,0.0, + half_thick_qua5 - half_thick_alu,0,"ONLY");  

  parbox[0] = half_thick_alu;
  parbox[1] = height5 / TMath::Cos(theta)/ 2.0;
  parbox[2] = half_thick_qua5;
  
  gMC->Gsvolu("V0E2","BOX",idtmed[3004],parbox,3);   // aluminum rectangular foil
  Float_t h5;
  h5 = TMath::Tan(theta) * (r0  + height1 + height2 + height3 + height4 + height5/2.0);  
  gMC->Gspos("V0E2",1,"V0R5",-h5 + half_thick_alu,0.0,0.0,idrotm[911],"ONLY"); 
  gMC->Gspos("V0E2",2,"V0R5",+h5 - half_thick_alu,0.0,0.0,idrotm[912],"ONLY"); 

  parbox[0] = TMath::Tan(theta) * (r0 + height1 + height2 + height3 + height4);
  parbox[1] = half_thick_alu;
  parbox[2] = half_thick_qua5;
  gMC->Gsvolu("V0E3","BOX",idtmed[3004],parbox,3);
  gMC->Gspos("V0E3",1,"V0R5",0.0, - (height5/2.0) + half_thick_alu ,0.0,0,"ONLY");
    

  parbox[0] = TMath::Tan(theta) * r5;
  parbox[1] = half_thick_alu;
  parbox[2] = half_thick_qua5;
  gMC->Gsvolu("V0E4","BOX",idtmed[3004],parbox,3);
  gMC->Gspos("V0E4",1,"V0R5",0.0, (height5/2.0) - half_thick_alu,0.0,0,"ONLY");


  Float_t  phi_rad ;  
  Float_t  phi_deg = 180./6.; 

 // Partie de droite : 
 
  for(Float_t  phi = 15.0; phi < 360.0; phi = phi + phi_deg)
      {
        phi_rad = phi*pi/180.;
      	AliMatrix(idrotm[902], 90.0, phi, 90.0, 90.0 +phi, 0.0 , 0.0);
        gMC->Gspos("V0R0",n_detec_R,"V0RI",-dist0*TMath::Sin(phi_rad),
	                  dist0*TMath::Cos(phi_rad),offset + half_thick_qua1,idrotm[902],"ONLY");
	n_detec_R++;
       }

  gMC->Gspos("V0RI",1,"ALIC",0.0,0.0,zdet,0,"ONLY");
 
  n_cells_R = (n_detec_R - 1) * 5;
  printf(" \n\n\n"); 
  printf("    Number of cells on Right side =   %d\n",  n_cells_R);    

// Partie de gauche :

  for(Float_t  phi = 15.0; phi < 360.0; phi = phi + phi_deg)
      {
        phi_rad = phi*pi/180.;
      	AliMatrix(idrotm[902], 90.0, phi, 90.0, 90.0 +phi, 0.0 , 0.0);
        gMC->Gspos("V0L0",n_detec_L,"V0LE",-dist0_left*TMath::Sin(phi_rad),
	                  dist0_left*TMath::Cos(phi_rad),offset_left + half_thick_qua1,idrotm[902],"ONLY");
        n_detec_L++;
       }


  gMC->Gspos("V0LE",1,"ALIC",0.0,0.0,-350.0-fThickness1/2.0,0,"ONLY");
 
  n_cells_L = (n_detec_L - 1) * 5;
  printf(" \n\n\n"); 
  printf("    Number of cells on Left side  =   %d\n\n",  n_cells_L);    

         
}
    
    
    
//_____________________________________________________________________________
void AliVZEROv0::BuildGeometry()
{
  
  // Builds simple ROOT TNode geometry for event display

  
  Int_t i;

  printf("\n");
  for(i=0;i<30;i++) printf("*");
  printf(" VZERO BuildGeometry ");
  for(i=0;i<30;i++) printf("*");
  printf("\n");
  
  TNode *Top; 

  TNode *V0Rnode, *V0Rnode0, *V0Rnode6 , *V0Rnode7, *V0Rnode8, *V0Rnode9;
  TNode *V0Rnode1, *V0Rnode2, *V0Rnode3, *V0Rnode4, *V0Rnode5;
  TNode *V0Lnode, *V0Lnode0;
  TNode *V0Lnode1, *V0Lnode2, *V0Lnode3, *V0Lnode4, *V0Lnode5;
   
  const int kColorVZERO  = kGreen;
 
  Top = gAlice->GetGeometry()->GetNode("alice");

  Float_t  height1, height2, height3, height4, height5; 
  Float_t  height;
  Float_t  theta;  
  Float_t  half_thick_alu;
  Float_t  half_thick_qua1,half_thick_qua2,half_thick_qua3;
  Float_t  half_thick_qua4,half_thick_qua5;
  Float_t  zdet;
  Float_t  r0, r5;
  Float_t  pi = TMath::Pi();
  Float_t  thick_alu;
  
//   height1           =     1.9;         
//   height2           =     3.7; 
//   height3           =     6.2;
//   height4           =    10.5;
//   height5           =    10.5;

  height1           =     2.0;           // height of cell 1, in cm
  height2           =     3.2;           // height of cell 2, in cm
  height3           =     4.9;           // height of cell 3, in cm
  height4           =     7.5;           // height of cell 4, in cm
  height5           =    12.0;           // height of cell 5, in cm
  
  theta             = pi/6.0/2.0;    
  half_thick_alu    = 0.0025;        
  thick_alu         = 2.0 * half_thick_alu; 
  half_thick_qua1   = fThickness1/2.0; 
  half_thick_qua2   = half_thick_qua1  - 0.25;
  half_thick_qua3   = half_thick_qua2  - 0.25;
  half_thick_qua4   = half_thick_qua3  - 0.25;
  half_thick_qua5   = half_thick_qua4  - 0.25;
  
  zdet              =    86.9 +fThickness/2.0;   
  r0                =    4.0;         
  height            =    height1 + height2 + height3 + height4 + height5;
  r5                =    r0 + height;

  Float_t   partube[3];

  partube[0] =  r0 - 0.2;
  partube[1] = (r5 + 1.0) / TMath::Cos(theta);
  partube[2] = fThickness/2.0;   
  
  TTUBE *V0RI = new TTUBE("V0RI", "V0RI", "void", partube[0], partube[1], partube[2]);
 		
  Top->cd();
  
  V0Rnode = new TNode("V0RI","V0RI",V0RI,0.0,0.0,+zdet,0);
  
  V0Rnode->SetLineColor(kBlue);
  fNodes->Add(V0Rnode);
  
  V0Rnode->SetVisibility(2);
     
 
// Rondelles de carbone (epaisseur 1 mm) de maintien des cellules ...
 

  Float_t   parbox[10];
 
  parbox[0] =    0.;
  parbox[1] =  360.;
  parbox[2] =   12.;
  parbox[3] =    2.;
  parbox[4] =  -0.1/2.0;
  parbox[5] =  r0;
  parbox[6] =  r5;     
  parbox[7] =  +0.1/2.0;
  parbox[8] =  r0;
  parbox[9] =  r5;  
  
  
  TPGON *V0CA = new TPGON("V0CA", "V0CA", "void",parbox[0], parbox[1],
			  static_cast<Int_t>(parbox[2]),
			  static_cast<Int_t>(parbox[3]));
			    
  V0CA->DefineSection( 0, parbox[4], parbox[5], parbox[6] );
  V0CA->DefineSection( 1, parbox[7], parbox[8], parbox[9] ); 
   
  V0Rnode->cd();
  V0Rnode6 = new TNode("V0CA", "V0CA",V0CA,0.0,0.0, fThickness/2.0-parbox[7],0);	 
  V0Rnode6->SetLineColor(kYellow);
  fNodes->Add(V0Rnode6); 
  V0Rnode->cd();
  V0Rnode7 = new TNode("V0CA", "V0CA",V0CA,0.0,0.0,-fThickness/2.0+parbox[7],0);	 
  V0Rnode7->SetLineColor(kYellow);
  fNodes->Add(V0Rnode7);
  
  parbox[4] =  -fThickness/2.0;
  parbox[5] =  r0 - 0.2;
  parbox[6] =  r0;     
  parbox[7] =  +fThickness/2.0;
  parbox[8] =  r0 - 0.2;
  parbox[9] =  r0; 

  TPGON  *V0IR = new TPGON("V0IR","V0IR","void",  parbox[0], parbox[1],
			  static_cast<Int_t>(parbox[2]),
			  static_cast<Int_t>(parbox[3]));
  V0IR->DefineSection( 0, parbox[4], parbox[5], parbox[6] );
  V0IR->DefineSection( 1, parbox[7], parbox[8], parbox[9] );
  
  V0Rnode->cd();
  V0Rnode8 = new TNode("V0IR", "V0IR",V0IR,0.0,0.0,0.0,0);
  V0Rnode8->SetLineColor(kYellow);
  fNodes->Add(V0Rnode8);
 
  parbox[4] =  -fThickness/2.0;
  parbox[5] =  r5;
  parbox[6] =  r5 + 1.0;     
  parbox[7] =  +fThickness/2.0;
  parbox[8] =  r5;
  parbox[9] =  r5 + 1.0; 

  TPGON  *V0ER = new TPGON("V0ER","V0ER","void",  parbox[0], parbox[1],
			  static_cast<Int_t>(parbox[2]),
			  static_cast<Int_t>(parbox[3]));

  V0ER->DefineSection( 0, parbox[4], parbox[5], parbox[6] );
  V0ER->DefineSection( 1, parbox[7], parbox[8], parbox[9] );
  
  V0Rnode->cd();
  V0Rnode9 = new TNode("V0ER", "V0ER",V0ER,0.0,0.0,0.0,0);
  V0Rnode9->SetLineColor(kYellow);
  fNodes->Add(V0Rnode9);

  Float_t   dist0; 
  Float_t   par[11];   
    
  dist0     =  r0 + height / 2.0;
  thick_alu =  2.0*half_thick_alu;
  
  par[0]    =  half_thick_qua1;
  par[1]    =  0.0;
  par[2]    =  0.0;
  par[3]    =  height / 2.0 ;
  par[4]    =  TMath::Tan(theta) * r0;
  par[5]    =  TMath::Tan(theta) * r5;
  par[6]    =  0.0;
  par[7]    =  height / 2.0 ;
  par[8]    =  TMath::Tan(theta) * r0;
  par[9]    =  TMath::Tan(theta) * r5;
  par[10]   =  0.0;
    
  TTRAP *V0R0 = new TTRAP("V0R0", "V0R0", "void", par[0], par[1], par[2], par[3],
  		par[4], par[5], par[6], par[7], par[8], par[9], par[10]);
		
  Float_t   dist1;
  Float_t   r1;
  Float_t   offset; 
    
  dist1     =  (- height + height1) /2.0; 
  r1        =  r0 + height1;
  offset    = - fThickness/2.0 + 0.1; 
   
  par[0]    =  half_thick_qua1 - thick_alu;
  par[3]    =  height1 / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r0 - thick_alu;
  par[5]    =  TMath::Tan(theta) * r1- thick_alu;
  par[7]    =  height1 / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r0 - thick_alu;
  par[9]    =  TMath::Tan(theta) * r1 - thick_alu;
  
  TTRAP *V0R1 = new TTRAP("V0R1", "V0R1", "void", par[0], par[1], par[2], par[3],
  		par[4], par[5], par[6], par[7], par[8], par[9], par[10]);
  		

  Float_t   dist2;
  Float_t   r2; 
    
  dist2     =    (- height + height2) /2.0 + height1;
  r2        =       r1 + height2; 
  
  par[0]    =  half_thick_qua2 - thick_alu;   
  par[3]    =  height2 / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r1 - thick_alu;
  par[5]    =  TMath::Tan(theta) * r2 - thick_alu;
  par[7]    =  height2 / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r1 - thick_alu;
  par[9]    =  TMath::Tan(theta) * r2 - thick_alu;


  TTRAP *V0R2 = new TTRAP("V0R2", "V0R2", "void", par[0], par[1], par[2], par[3],
  		par[4], par[5], par[6], par[7], par[8], par[9], par[10]);
 

  Float_t   dist3;
  Float_t   r3; 
    
  dist3     =    (- height + height3) /2.0 + height1 + height2;
  r3        =       r2 + height3; 
   
  par[0]    =  half_thick_qua3 - thick_alu;   
  par[3]    =  height3 / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r2 - thick_alu;
  par[5]    =  TMath::Tan(theta) * r3 - thick_alu;
  par[7]    =  height3 / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r2 - thick_alu;
  par[9]    =  TMath::Tan(theta) * r3 - thick_alu;


  TTRAP *V0R3 = new TTRAP("V0R3", "V0R3", "void", par[0], par[1], par[2], par[3],
  		par[4], par[5], par[6], par[7], par[8], par[9], par[10]);
 

  Float_t   dist4;
  Float_t   r4; 
    
  dist4     =    (- height + height4) /2.0 + height1 + height2 + height3;
  r4        =       r3 + height4; 
   
  par[0]    =  half_thick_qua4 - thick_alu;   
  par[3]    =  height4 / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r3 - thick_alu;
  par[5]    =  TMath::Tan(theta) * r4 - thick_alu;
  par[7]    =  height4 / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r3 - thick_alu;
  par[9]    =  TMath::Tan(theta) * r4 - thick_alu;


  TTRAP *V0R4 = new TTRAP("V0R4", "V0R4", "void", par[0], par[1], par[2], par[3],
  		par[4], par[5], par[6], par[7], par[8], par[9], par[10]);


  Float_t   dist5;
    
  dist5     =    (- height + height5) /2.0 + height1 + height2 + height3 + height4;
  
  par[0]    =  half_thick_qua5 - thick_alu;      
  par[3]    =  height5 / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r4 - thick_alu;
  par[5]    =  TMath::Tan(theta) * r5 - thick_alu;
  par[7]    =  height5 / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r4 - thick_alu;
  par[9]    =  TMath::Tan(theta) * r5 - thick_alu;


  TTRAP *V0R5 = new TTRAP("V0R5", "V0R5", "void", par[0], par[1], par[2], par[3],
  		par[4], par[5], par[6], par[7], par[8], par[9], par[10]);

		
  Float_t  phi;
  Float_t  phi_deg= 180./6.;
  Float_t  phi_rad;
  Float_t  xdet,ydet;
  Int_t    n_detec_R = 1; 

  char     NameNode[12];  

 
  for (phi = 15.0; phi < 360.0; phi = phi + phi_deg) 
  {
     
    TRotMatrix* mat920 = new TRotMatrix("rot920","rot920", 90.0, +phi, 90., 90.+phi, 0.0, 0.0 );        
    
    phi_rad = phi*pi/180.;
    xdet = dist0*TMath::Sin(phi_rad);
    ydet = dist0*TMath::Cos(phi_rad);
    
  
    sprintf(NameNode,"SUBDER%d",n_detec_R);
    
    V0Rnode->cd();
    V0Rnode0 = new TNode(NameNode,NameNode,V0R0,-xdet,ydet, offset + half_thick_qua1,mat920);	 
    V0Rnode0->SetLineColor(kColorVZERO);
    fNodes->Add(V0Rnode0);
    n_detec_R++;
    
    sprintf(NameNode,"SUBDER%d",n_detec_R);
    V0Rnode0->cd();    
    V0Rnode1 = new TNode(NameNode,NameNode,V0R1,0.0,dist1, 0.0,0);	 
    V0Rnode1->SetLineColor(kColorVZERO);
    fNodes->Add(V0Rnode1);
    n_detec_R++;
    
    sprintf(NameNode,"SUBDER%d",n_detec_R);
    V0Rnode0->cd();    
    V0Rnode2 = new TNode(NameNode,NameNode,V0R2,0.0,dist2, - half_thick_qua1 + half_thick_qua2,0);	 
    V0Rnode2->SetLineColor(kColorVZERO);
    fNodes->Add(V0Rnode2);
    n_detec_R++;


    sprintf(NameNode,"SUBDER%d",n_detec_R);
    V0Rnode0->cd();    
    V0Rnode3 = new TNode(NameNode,NameNode,V0R3,0.0,dist3, - half_thick_qua1 + half_thick_qua3,0);	 
    V0Rnode3->SetLineColor(kColorVZERO);
    fNodes->Add(V0Rnode3);
    n_detec_R++;

    sprintf(NameNode,"SUBDER%d",n_detec_R);
    V0Rnode0->cd();    
    V0Rnode4 = new TNode(NameNode,NameNode,V0R4,0.0,dist4, - half_thick_qua1 + half_thick_qua4,0);	 
    V0Rnode4->SetLineColor(kColorVZERO);
    fNodes->Add(V0Rnode4);
    n_detec_R++;
     
    sprintf(NameNode,"SUBDER%d",n_detec_R);
    V0Rnode0->cd();    
    V0Rnode5 = new TNode(NameNode,NameNode,V0R5,0.0,dist5, - half_thick_qua1 + half_thick_qua5,0);	 
    V0Rnode5->SetLineColor(kColorVZERO);
    fNodes->Add(V0Rnode5);
    n_detec_R++;
       
    V0Rnode0->SetVisibility(2);
    
  }  


// Left side of VZERO :
  
  
  Float_t   dist0_left;   
  Float_t   r0_left      =   4.3;   
  Float_t   height1_left =   2.6; 
  Float_t   height2_left =   4.1;
  Float_t   height3_left =   6.4;
  Float_t   height4_left =  10.2;
  Float_t   height5_left =  16.9;
  Float_t   height_left  = height1_left + height2_left + height3_left 
                                        + height4_left + height5_left;
  Float_t   r5_left      = r0_left  + height_left; 

  partube[0] =  r0_left;
  partube[1] = (r5_left) / TMath::Cos(theta);
  partube[2] = fThickness1/2.0; 
  
  TTUBE *V0LE = new TTUBE("V0LE", "V0LE", "void", partube[0], partube[1], partube[2]);
 		
  Top->cd();
  
  V0Lnode = new TNode("V0LE","V0LE",V0LE,0.0,0.0,-350.0-fThickness1/2.0,0);
  
  V0Lnode->SetLineColor(kBlue);
  fNodes->Add(V0Lnode);
  
  V0Lnode->SetVisibility(2);
  
  dist0_left  =  r0_left + height_left / 2.0;
  thick_alu   =  2.0*half_thick_alu;
  
  par[0]      =  half_thick_qua1;
  par[1]      =  0.0;
  par[2]      =  0.0;
  par[3]      =  height_left / 2.0 ;
  par[4]      =  TMath::Tan(theta) * r0_left;
  par[5]      =  TMath::Tan(theta) * r5_left;
  par[6]      =  0.0;
  par[7]      =  height_left / 2.0 ;
  par[8]      =  TMath::Tan(theta) * r0_left;
  par[9]      =  TMath::Tan(theta) * r5_left;
  par[10]     =  0.0;
  
  TTRAP *V0L0 = new TTRAP("V0L0", "V0L0", "void", par[0], par[1], par[2], par[3],
  		par[4], par[5], par[6], par[7], par[8], par[9], par[10]);
		
  
  Float_t   dist1_left;
  Float_t   r1_left;
  Float_t   offset_left;
     
  dist1_left     =  (- height_left + height1_left) /2.0; 
  r1_left        =  r0_left + height1_left;
  offset_left    = - fThickness1/2.0 + 0.1; 
   
  par[0]    =  half_thick_qua1 - thick_alu;
  par[3]    =  height1_left / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r0_left - thick_alu;
  par[5]    =  TMath::Tan(theta) * r1_left- thick_alu;
  par[7]    =  height1_left / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r0_left - thick_alu;
  par[9]    =  TMath::Tan(theta) * r1_left - thick_alu;

  TTRAP *V0L1 = new TTRAP("V0L1", "V0L1", "void", par[0], par[1], par[2], par[3],
  		par[4], par[5], par[6], par[7], par[8], par[9], par[10]);
 
  Float_t   dist2_left;
  Float_t   r2_left; 
    
  dist2_left     =    (- height_left + height2_left) /2.0 + height1_left;
  r2_left        =       r1_left + height2_left; 
  
  par[0]    =  half_thick_qua1 - thick_alu;   
  par[3]    =  height2_left / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r1_left - thick_alu;
  par[5]    =  TMath::Tan(theta) * r2_left - thick_alu;
  par[7]    =  height2_left / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r1_left - thick_alu;
  par[9]    =  TMath::Tan(theta) * r2_left - thick_alu;

  TTRAP *V0L2 = new TTRAP("V0L2", "V0L2", "void", par[0], par[1], par[2], par[3],
  		par[4], par[5], par[6], par[7], par[8], par[9], par[10]);
 

  
  Float_t   dist3_left;
  Float_t   r3_left;
     
  dist3_left     =    (- height_left + height3_left) /2.0 + height1_left + height2_left;
  r3_left        =       r2_left + height3_left; 
   
  par[0]    =  half_thick_qua1 - thick_alu;   
  par[3]    =  height3_left / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r2_left - thick_alu;
  par[5]    =  TMath::Tan(theta) * r3_left - thick_alu;
  par[7]    =  height3_left / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r2_left - thick_alu;
  par[9]    =  TMath::Tan(theta) * r3_left - thick_alu;
  
  TTRAP *V0L3 = new TTRAP("V0L3", "V0L3", "void", par[0], par[1], par[2], par[3],
  		par[4], par[5], par[6], par[7], par[8], par[9], par[10]);
 
   
  Float_t   dist4_left;
  Float_t   r4_left;
       
  dist4_left     =    (- height_left + height4_left) /2.0 + height1_left 
                                     + height2_left + height3_left;
  r4_left        =       r3_left + height4_left; 
   
  par[0]    =  half_thick_qua1 - thick_alu;   
  par[3]    =  height4_left / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r3_left - thick_alu;
  par[5]    =  TMath::Tan(theta) * r4_left - thick_alu;
  par[7]    =  height4_left / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r3_left - thick_alu;
  par[9]    =  TMath::Tan(theta) * r4_left - thick_alu;
  
  TTRAP *V0L4 = new TTRAP("V0L4", "V0L4", "void", par[0], par[1], par[2], par[3],
  		par[4], par[5], par[6], par[7], par[8], par[9], par[10]);
		

  Float_t   dist5_left;

         
  dist5_left     =    (- height_left + height5_left) /2.0 + height1_left 
                                     + height2_left + height3_left + height4_left;

   
  par[0]    =  half_thick_qua1 - thick_alu;   
  par[3]    =  height5_left / 2.0 - thick_alu;
  par[4]    =  TMath::Tan(theta) * r4_left - thick_alu;
  par[5]    =  TMath::Tan(theta) * r5_left - thick_alu;
  par[7]    =  height5_left / 2.0 - thick_alu;
  par[8]    =  TMath::Tan(theta) * r4_left - thick_alu;
  par[9]    =  TMath::Tan(theta) * r5_left - thick_alu;
  
  TTRAP *V0L5 = new TTRAP("V0L5", "V0L5", "void", par[0], par[1], par[2], par[3],
  		par[4], par[5], par[6], par[7], par[8], par[9], par[10]);


  Int_t    n_detec_L = 1;
 
  for (phi = 15.0; phi < 360.0; phi = phi + phi_deg) 
  {
     
    TRotMatrix* mat920 = new TRotMatrix("rot920","rot920", 90.0, +phi, 90., 90.+phi, 0.0, 0.0 );        
    
    phi_rad = phi*pi/180.;
    xdet = dist0_left*TMath::Sin(phi_rad);
    ydet = dist0_left*TMath::Cos(phi_rad);
    
  
    sprintf(NameNode,"SUBDEL%d",n_detec_L);
    
    V0Lnode->cd();
    V0Lnode0 = new TNode(NameNode,NameNode,V0L0,-xdet,ydet, offset_left + half_thick_qua1,mat920);	 
    V0Lnode0->SetLineColor(kColorVZERO);
    fNodes->Add(V0Lnode0);
    n_detec_L++;
    
    sprintf(NameNode,"SUBDEL%d",n_detec_L);
    V0Lnode0->cd();    
    V0Lnode1 = new TNode(NameNode,NameNode,V0L1,0.0,dist1_left, 0.0,0);	 
    V0Lnode1->SetLineColor(kColorVZERO);
    fNodes->Add(V0Lnode1);
    n_detec_L++;
    
    sprintf(NameNode,"SUBDEL%d",n_detec_L);
    V0Lnode0->cd();    
    V0Lnode2 = new TNode(NameNode,NameNode,V0L2,0.0,dist2_left, 0.0,0);	 
    V0Lnode2->SetLineColor(kColorVZERO);
    fNodes->Add(V0Lnode2);
    n_detec_L++;


    sprintf(NameNode,"SUBDEL%d",n_detec_L);
    V0Lnode0->cd();    
    V0Lnode3 = new TNode(NameNode,NameNode,V0L3,0.0,dist3_left, 0.0,0);	 
    V0Lnode3->SetLineColor(kColorVZERO);
    fNodes->Add(V0Lnode3);
    n_detec_L++;

    sprintf(NameNode,"SUBDEL%d",n_detec_L);
    V0Lnode0->cd();    
    V0Lnode4 = new TNode(NameNode,NameNode,V0L4,0.0,dist4_left, 0.0,0);	 
    V0Lnode4->SetLineColor(kColorVZERO);
    fNodes->Add(V0Lnode4);
    n_detec_L++;
     
    sprintf(NameNode,"SUBDEL%d",n_detec_L);
    V0Lnode0->cd();    
    V0Lnode5 = new TNode(NameNode,NameNode,V0L5,0.0,dist5_left, 0.0,0);	 
    V0Lnode5->SetLineColor(kColorVZERO);
    fNodes->Add(V0Lnode5);
    n_detec_L++;
       
    V0Lnode0->SetVisibility(2);
    
  }    

        
}  
    

//------------------------------------------------------------------------
void AliVZEROv0::CreateMaterials()
{
    Int_t i;

    printf("\n");
    for(i=0;i<30;i++) printf("*");
    printf(" VZERO create materials ");
    for(i=0;i<30;i++) printf("*");
    printf("\n");
    

    Float_t ppckov[14] = { 5.5e-9, 5.7e-9, 5.9e-9, 6.1e-9, 6.3e-9, 6.5e-9, 6.7e-9, 
                           6.9e-9, 7.1e-9, 7.3e-9, 7.5e-9, 7.7e-9, 7.9e-9, 8.1e-9 };

           
    Float_t ppckov_alu[14] = { 5.5e-9, 5.7e-9, 5.9e-9, 6.1e-9, 6.3e-9, 6.5e-9, 6.7e-9, 
                               6.9e-9, 7.1e-9, 7.3e-9, 7.5e-9, 7.7e-9, 7.9e-9, 8.1e-9 };
			   
    Float_t rindex_quarz[14] = { 1.52398,  1.53090, 1.53835, 1.54641, 1.55513, 1.56458, 
                                 1.57488,  1.58611, 1.59842, 1.61197, 1.62696, 1.64362, 
                                 1.662295, 1.68337 };
				 
    Float_t absco_quarz[14] = { 105.8,  45.656, 35.665, 28.598, 25.007, 21.04, 17.525, 
                                14.177, 9.282, 4.0925, 1.149, 0.3627, 0.1497, 0.05 }; 	
				  							
    Float_t effic_all[14]   = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    
        
    Float_t rindex_alu[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. }; 
    
    
    Float_t absco_alu[14]  = { 1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,
			       1e-4,1e-4,1e-4,1e-4 };
    Float_t effic_alu[14]  = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };


    Int_t *idtmed = fIdtmed->GetArray()-2999;
    
    
//  Parameters related to Quarz (SiO2) :
 
    Float_t aqua[2], zqua[2], densqua, wmatqua[2];
    Int_t nlmatqua;
    
    aqua[0]    = 28.09;
    aqua[1]    = 16.;
    zqua[0]    = 14.;
    zqua[1]    = 8.;
    densqua    = 2.64;
    nlmatqua   = -2;
    wmatqua[0] = 1.;
    wmatqua[1] = 2.;

// Parameters  related to aluminum sheets :
    
    Float_t  aal   = 26.98;
    Float_t  zal   = 13.00; 
    Float_t  densal=   2.7; 
    Float_t  radlal=   8.9;
       
// Parameters  related to scintillator CH :
    
    Float_t ascin[2] = {1.01,12.01};
    Float_t zscin[2] = {1,6};
    Float_t wscin[2] = {1,1};
    Float_t denscin  = 1.03;
    
//  Definition of materials :
       
    AliMaterial( 1, "AIR A$", 14.61, 7.3, .001205, 30420., 67500, 0, 0);
    AliMaterial(11, "AIR I$", 14.61, 7.3, .001205, 30420., 67500, 0, 0);
    AliMaterial( 2, "CARBON$"  , 12.01, 6.0, 2.265, 18.8, 49.9, 0, 0);
    AliMixture(  3, "QUA", aqua, zqua, densqua, nlmatqua, wmatqua);
    AliMaterial( 4, "ALUMINIUM1$", 26.98, 13., 2.7, 8.9, 37.2, 0, 0);
    AliMaterial( 5, "ALUMINIUM2$", aal, zal, densal, radlal, 0, 0, 0);
    
    
    AliMixture( 6, "Scintillator$",ascin,zscin,denscin,-2,wscin);
    
     
    Int_t   ISXFLD = gAlice->Field()->Integ();
    Float_t SXMGMX = gAlice->Field()->Max();
    
    Float_t tmaxfd, stemax, deemax, epsil, stmin;
    
    tmaxfd = 10.;
    stemax = 0.1;
    deemax = 0.1;     
    epsil  = 0.001;
    stmin  = 0.001;
              
    printf(" \n");
    printf(" StepQua,    StepAlu    = %f %f \n",fMaxStepQua,fMaxStepAlu);
    printf(" DeStepQua,  DeStepAlu  = %f %f \n",fMaxDestepQua,fMaxDestepAlu);
    printf(" \n");    


//  Active Air :    
    AliMedium(1, "ACTIVE AIR$", 1, 1, ISXFLD, SXMGMX,
              10.0, 1.0, 0.1, 0.1, 10.0, 0, 0) ;

//  Inactive air : 
  
    AliMedium(11, "INACTIVE AIR$", 11, 0, ISXFLD, SXMGMX,
              10.0, 1.0, 0.1, 0.1, 10.0, 0, 0) ;
    
    AliMedium(2, "CARBON$ ", 2,  1, ISXFLD, SXMGMX,
              tmaxfd, stemax, deemax, epsil, stmin, 0, 0);   

    AliMedium(3, "QUARZ$", 3, 1, ISXFLD, SXMGMX,
              tmaxfd, fMaxStepQua, fMaxDestepQua, epsil, stmin, 0, 0);
    
    AliMedium(4,"ALUMINUM1$",4, 1, ISXFLD, SXMGMX,
              tmaxfd, fMaxStepAlu, fMaxDestepAlu, epsil, stmin, 0, 0);
	      

    AliMedium(5,"ALUMINUM2$",5, 1, ISXFLD, SXMGMX,
              tmaxfd, fMaxStepAlu, fMaxDestepAlu, epsil, stmin, 0, 0);    

    AliMedium(6,"SCINTILLATOR$",6, 1, ISXFLD, SXMGMX, 10., .01, 1., .003, .003, 0, 0);

    gMC->Gstpar(idtmed[3000], "LOSS", 1.);  //  [3000] = air ACTIF  [3010] = air INACTIF
    gMC->Gstpar(idtmed[3000], "HADR", 1.);
    gMC->Gstpar(idtmed[3000], "DCAY", 1.);
    gMC->Gstpar(idtmed[3000], "DRAY", 1.);
    
    gMC->Gstpar(idtmed[3001], "LOSS", 1.);  //  [3001] = carbon
    gMC->Gstpar(idtmed[3001], "HADR", 1.);
    gMC->Gstpar(idtmed[3001], "DCAY", 1.);
    gMC->Gstpar(idtmed[3001], "DRAY", 1.);

    gMC->Gstpar(idtmed[3002], "LOSS", 1.);  //  [3002] = quartz
    gMC->Gstpar(idtmed[3002], "HADR", 1.);
    gMC->Gstpar(idtmed[3002], "DCAY", 1.);
    gMC->Gstpar(idtmed[3002], "DRAY", 1.);  
    gMC->Gstpar(idtmed[3002], "CUTGAM",0.5E-4) ; 
    gMC->Gstpar(idtmed[3002], "CUTELE",1.0E-4) ;
    
    gMC->Gstpar(idtmed[3003], "LOSS", 1.);  //  [3003] = normal aluminum
    gMC->Gstpar(idtmed[3003], "HADR", 1.);
    gMC->Gstpar(idtmed[3003], "DCAY", 1.);
    gMC->Gstpar(idtmed[3003], "DRAY", 1.);
    
    gMC->Gstpar(idtmed[3004], "LOSS", 1.);  //  [3004] = reflecting aluminum
    gMC->Gstpar(idtmed[3004], "HADR", 1.);
    gMC->Gstpar(idtmed[3004], "DCAY", 1.);
    gMC->Gstpar(idtmed[3004], "DRAY", 1.);
    gMC->Gstpar(idtmed[3004], "CUTGAM",0.5E-4) ; 
    gMC->Gstpar(idtmed[3004], "CUTELE",1.0E-4) ;
    
    gMC->Gstpar(idtmed[3005], "LOSS", 1.);  //  [3005] = scintillator
    gMC->Gstpar(idtmed[3005], "HADR", 1.);
    gMC->Gstpar(idtmed[3005], "DCAY", 1.);
    gMC->Gstpar(idtmed[3005], "DRAY", 1.);    
    
    gMC->SetCerenkov(idtmed[3002], 14, ppckov, absco_quarz, effic_all,rindex_quarz);    
    gMC->SetCerenkov(idtmed[3004], 14, ppckov_alu, absco_alu, effic_alu, rindex_alu);

    
}
//---------------------------------------------------------------------
void AliVZEROv0::DrawModule()
{

//  Drawing is done in DrawVZERO.C

   Int_t i;

   printf("\n");
   for(i=0;i<30;i++) printf("*");
   printf(" VZERO DrawModule ");
   for(i=0;i<30;i++) printf("*");
   printf("\n");


}

//-------------------------------------------------------------------
void AliVZEROv0::Init()
{
// Initialises version 0 of the VZERO Detector
// Just prints an information message
  
   printf(" VZERO version %d initialized \n",IsVersion());
   
//   gMC->SetMaxStep(fMaxStepAlu);
//   gMC->SetMaxStep(fMaxStepQua);
   
//   AliVZERO::Init();
  
}

//-------------------------------------------------------------------

void AliVZEROv0::StepManager()
{
  
//   (Very)Minimal version of StepManager 

    
     Int_t  copy;
     static Int_t vol[4];
     static Float_t hits[15];
     
     TLorentzVector pos;
     
     TLorentzVector mom;
     Float_t        theta;
     Float_t        phi;
     Float_t        kRaddeg = 180/TMath::Pi();
     Float_t        RingNumber;

     Int_t ipart;
     
          
//     TGeant3 *geant3 = (TGeant3*) gMC;     
//     Int_t  Nphot = geant3->Gckin2()->ngphot;
     

//   Only charged tracks :
     
     if ( !gMC->TrackCharge() || !gMC->IsTrackAlive() ) return; 


     vol[0]    = gMC->CurrentVolOffID(1, vol[1]);
     vol[2]    = gMC->CurrentVolID(copy);
     vol[3]    = copy;

     if      ( gMC->CurrentVolID(copy) == gMC->VolId("V0R1") ||
               gMC->CurrentVolID(copy) == gMC->VolId("V0L1") )
	       RingNumber = 1.0;
     else if ( gMC->CurrentVolID(copy) == gMC->VolId("V0R2") ||
               gMC->CurrentVolID(copy) == gMC->VolId("V0L2") ) 
	       RingNumber = 2.0;  
     else if ( gMC->CurrentVolID(copy) == gMC->VolId("V0R3") ||
               gMC->CurrentVolID(copy) == gMC->VolId("V0L3") )
	       RingNumber = 3.0;
     else if ( gMC->CurrentVolID(copy) == gMC->VolId("V0R4") ||
               gMC->CurrentVolID(copy) == gMC->VolId("V0L4") ) 	 
	       RingNumber = 4.0; 
     else if ( gMC->CurrentVolID(copy) == gMC->VolId("V0R5") ||
               gMC->CurrentVolID(copy) == gMC->VolId("V0L5") )	  
               RingNumber = 5.0; 
     else
     	       RingNumber = 0.0;

     if (gMC->IsTrackEntering() && RingNumber > 0.5) {
       
         gMC->TrackPosition(pos);
     
         gMC->TrackMomentum(mom);      
         Double_t tc   = mom[0]*mom[0]+mom[1]*mom[1];
         Double_t Pt   = TMath::Sqrt(tc);
	 Double_t Pmom = TMath::Sqrt(tc+mom[2]*mom[2]);

         theta   = Float_t(TMath::ATan2(Pt,Double_t(mom[2])))*kRaddeg;
         phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
     
         ipart  = gMC->TrackPid();

         hits[0]  = pos[0];
         hits[1]  = pos[1];
         hits[2]  = pos[2];	 	 
	 hits[3]  =  ipart; 
	 
//         Float_t ttime = gMC->TrackTime();
//         hits[4] = ttime*1e9;

	 hits[4]  = gMC->TrackTime();
         hits[5]  = gMC->TrackCharge();
	 hits[6]  = theta;
	 hits[7]  = phi;
	 hits[8]  = RingNumber;
	 
	 hits[9]  = Pt;
	 hits[10] = Pmom;
	 hits[11] = mom[0];
	 hits[12] = mom[1];
	 hits[13] = mom[2];
	 

         AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);
	 
	 }
     
}

//_____________________________________________________________________________
void AliVZEROv0::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a VZERO hit
  //

  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliVZEROhit(fIshunt,track,vol,hits);
}

//---------------------------------------------------------------------
void AliVZEROv0::AddDigits(Int_t *tracks, Int_t* digits) 
{

   TClonesArray  &ldigits = *fDigits;
   new(ldigits[fNdigits++]) AliVZEROdigit(tracks, digits);
}

//---------------------------------------------------------------------
void AliVZEROv0::MakeBranch(Option_t *option)
{
  
  // Creates new branches in the current Root Tree
  
  
  char branchname[10];
  sprintf(branchname,"%s",GetName());
  printf(" fBufferSize = %d \n",fBufferSize);
  
  const char *H = strstr(option,"H");
  
  if (fHits   && TreeH() && H) {
    TreeH()->Branch(branchname,&fHits, fBufferSize);
    printf("* AliDetector::MakeBranch * Making Branch %s for hits\n",branchname);
  }     

  const char *D = strstr(option,"D");
  //
  if (fDigits   && gAlice->TreeD() && D) {
    gAlice->TreeD()->Branch(branchname,&fDigits, fBufferSize);
    printf("* AliDetector::MakeBranch * Making Branch %s for digits\n",branchname);
  }  
   
}
