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



//////////////////////////////////////////////////////////////////////
//                                                                  //
//  (V-zero) detector  version 0  as designed by the Lyon group     //
//   All comments should be sent to Brigitte CHEYNIS :              //
//                                  b.cheynis@ipnl.in2p3.fr         // 
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TGeometry.h>
#include <TTRD2.h>
#include <TCONE.h>
#include <TPGON.h>
#include <TPCON.h>
#include <TSPHE.h>
#include <TTRAP.h>
#include <TBRIK.h>
#include <TBox.h>

#include <TShape.h>
#include <TNode.h>
#include <TClonesArray.h>
#include <TH1.h>
#include <string.h>
#include <iostream.h>

#include "AliVZEROv0.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliVZEROhit.h"
#include "AliVZEROdigit.h"
#include <iostream.h>
#include <fstream.h>

#include <TGeant3.h>
#include <stdlib.h>
#include "TObjectTable.h"

#include "AliConst.h"
#include "ABSOSHILConst.h"
#include "ABSOConst.h"

ClassImp(AliVZEROv0)

//--------------------------------------------------------------------
AliVZEROv0:: AliVZEROv0():AliVZERO()
{
  fRootFile = 0;
  fhMultiplicity = 0;
  fhGEANTcode = 0;
  fhCerenkov = 0;
  fhToF = 0;
}
//--------------------------------------------------------------------
AliVZEROv0::AliVZEROv0(const char *name, const char *title):
 AliVZERO(name,title)
{

// Standard constructor for V-zeroR Detector (right part)  version 0


  Int_t i;

  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" create VZERO object");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  
}
//-------------------------------------------------------------------------
void AliVZEROv0::CreateGeometry()
{

// Creates the Geant geometry of the V-zeroR Detector (right part) version 0

  
  Int_t i;
  
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" VZERO Create Geometry ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
    
  
  Int_t    *idtmed = fIdtmed->GetArray()-2999;

  Int_t    n_detec = 1;
  Int_t    n_cells = 1;

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
  
  height1           =     2.3;           // height of cell 1, in cm
  height2           =     3.7;           // height of cell 2, in cm
  height3           =     6.2;           // height of cell 3, in cm
  height4           =    10.5;           // height of cell 4, in cm
  height5           =    10.5;           // height of cell 5, in cm
  
  theta             = pi/12.0/2.0;       // half angular opening = 7.5 degrees
  half_thick_alu    = 0.0025;            // half thickness of aluminum foil, in cm
  thick_alu         = 2.0 * half_thick_alu; 
  half_thick_qua1   = fThickness1/2.0;   // half thickness of WRAPPED quartz cell (inner ring)
  half_thick_qua2   = half_thick_qua1  - 0.25;
  half_thick_qua3   = half_thick_qua2  - 0.25;
  half_thick_qua4   = half_thick_qua3  - 0.25;
  half_thick_qua5   = half_thick_qua4  - 0.25;
  
  zdet              =    86.9 +fThickness/2.0;  // distance to vertex (along Z axis)
  r0                =    3.4;            // closest distance to center of the beam pipe
  height            =    height1 + height2 + height3 + height4 + height5;
  r5                =    r0 + height;


// Creation of mother volume V0RI :

  Float_t   partube[3];
  
  partube[0] = r0 - 0.2;
  partube[1] = (r5 + 1.0) / TMath::Cos(theta);
  partube[2] = fThickness/2.0; 
  
    
  gMC->Gsvolu("V0RI","TUBE",idtmed[3010],partube,3);
  
// Creation of  carbon lids (1 mm thick) to keep V0RI box shut...

  Float_t   parbox[10];
  
  parbox[0] =    0.;
  parbox[1] =  360.;
  parbox[2] =    24;
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
  Float_t  theta_deg = 180./12./2.0; 
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
  Float_t  phi_deg = 180./12.; 

  
  for(Float_t  phi = 7.5; phi < 360.0; phi = phi + phi_deg)
      {
        phi_rad = phi*pi/180.;
      	AliMatrix(idrotm[902], 90.0, phi, 90.0, 90.0 +phi, 0.0 , 0.0);
        gMC->Gspos("V0R0",n_detec,"V0RI",-dist0*TMath::Sin(phi_rad),
	                  dist0*TMath::Cos(phi_rad),offset + half_thick_qua1,idrotm[902],"ONLY");
	n_detec++;
       }


  gMC->Gspos("V0RI",1,"alic",0.0,0.0,zdet,0,"ONLY");
 
  n_cells = (n_detec - 1) * 5;
  printf(" \n\n\n"); 
  printf(" Number of cells =   %d\n\n",  n_cells);    
         
}
    
    
//_____________________________________________________________________________
void AliVZEROv0::BuildGeometry()
{
  
  // Builds simple ROOT TNode geometry for event display
  
  
  Int_t i;

  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" VZERO BuildGeometry ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  
  TNode *Top, *Node0, *Node1, *Node2;
  TNode *Node3, *Node4 , *Node5, *Node6 , *Node7;
  TNode *V0Rnode, *V0Rnode0, *V0Rnode6 , *V0Rnode7, *V0Rnode8, *V0Rnode9;
  TNode *V0Rnode1, *V0Rnode2, *V0Rnode3, *V0Rnode4, *V0Rnode5;
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
  
  height1           =     2.3;         
  height2           =     3.7; 
  height3           =     6.2;
  height4           =    10.5;
  height5           =    10.5;
  
  theta             = pi/12.0/2.0;    
  half_thick_alu    = 0.0025;        
  thick_alu         = 2.0 * half_thick_alu; 
  half_thick_qua1   = fThickness1/2.0; 
  half_thick_qua2   = half_thick_qua1  - 0.25;
  half_thick_qua3   = half_thick_qua2  - 0.25;
  half_thick_qua4   = half_thick_qua3  - 0.25;
  half_thick_qua5   = half_thick_qua4  - 0.25;
  
  zdet              =    86.9 +fThickness/2.0;   
  r0                =    3.4;         
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
  parbox[2] =    24;
  parbox[3] =     2;
  parbox[4] =  -0.1/2.0;
  parbox[5] =  r0;
  parbox[6] =  r5;     
  parbox[7] =  +0.1/2.0;
  parbox[8] =  r0;
  parbox[9] =  r5;  
  
  
  TPGON *V0CA = new TPGON("V0CA", "V0CA", "void",parbox[0], parbox[1],
                            parbox[2],parbox[3]);
			    
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
                            parbox[2],parbox[3]);
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
                            parbox[2],parbox[3]);
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
  Float_t  phi_deg= 180./12.;
  Float_t  phi_rad;
  Float_t  xdet,ydet;
  Int_t    n_detec = 1; 
  char     NameNode[12];  

 
  for (phi = 7.5; phi < 360.0; phi = phi + phi_deg) 
  {
     
    TRotMatrix* mat920 = new TRotMatrix("rot920","rot920", 90.0, +phi, 90., 90.+phi, 0.0, 0.0 );        
    
    phi_rad = phi*pi/180.;
    xdet = dist0*TMath::Sin(phi_rad);
    ydet = dist0*TMath::Cos(phi_rad);
    
  
    sprintf(NameNode,"SUBDET%d",n_detec);
    
    V0Rnode->cd();
    V0Rnode0 = new TNode(NameNode,NameNode,V0R0,-xdet,ydet, offset + half_thick_qua1,mat920);	 
    V0Rnode0->SetLineColor(kColorVZERO);
    fNodes->Add(V0Rnode0);
    n_detec++;
    
    sprintf(NameNode,"SUBDET%d",n_detec);
    V0Rnode0->cd();    
    V0Rnode1 = new TNode(NameNode,NameNode,V0R1,0.0,dist1, 0.0,0);	 
    V0Rnode1->SetLineColor(kColorVZERO);
    fNodes->Add(V0Rnode1);
    n_detec++;
    
    sprintf(NameNode,"SUBDET%d",n_detec);
    V0Rnode0->cd();    
    V0Rnode2 = new TNode(NameNode,NameNode,V0R2,0.0,dist2, - half_thick_qua1 + half_thick_qua2,0);	 
    V0Rnode2->SetLineColor(kColorVZERO);
    fNodes->Add(V0Rnode2);
    n_detec++;


    sprintf(NameNode,"SUBDET%d",n_detec);
    V0Rnode0->cd();    
    V0Rnode3 = new TNode(NameNode,NameNode,V0R3,0.0,dist3, - half_thick_qua1 + half_thick_qua3,0);	 
    V0Rnode3->SetLineColor(kColorVZERO);
    fNodes->Add(V0Rnode3);
    n_detec++;

    sprintf(NameNode,"SUBDET%d",n_detec);
    V0Rnode0->cd();    
    V0Rnode4 = new TNode(NameNode,NameNode,V0R4,0.0,dist4, - half_thick_qua1 + half_thick_qua4,0);	 
    V0Rnode4->SetLineColor(kColorVZERO);
    fNodes->Add(V0Rnode4);
    n_detec++;
     
    sprintf(NameNode,"SUBDET%d",n_detec);
    V0Rnode0->cd();    
    V0Rnode5 = new TNode(NameNode,NameNode,V0R5,0.0,dist5, - half_thick_qua1 + half_thick_qua5,0);	 
    V0Rnode5->SetLineColor(kColorVZERO);
    fNodes->Add(V0Rnode5);
    n_detec++;
       
    V0Rnode0->SetVisibility(2);
    
  }  
  

//  Here I add the enveloppe of the  FRONTAL ABSORBER defined by Andreas Morsch 
//  in   AliABSOv0::CreateGeometry()  :

    Float_t parm[24];
    Float_t dz;
    
    Top->cd();
  
    parm[0]  = 0.;
    parm[1]  = 360.;
    parm[2]  = 7.;

    parm[3]  = -(zRear-zAbsStart)/2.;
    parm[4]  = rAbs;
    parm[5]  = zAbsStart * TMath::Tan(theta1);

    parm[6]  = parm[3]+(zNose-zAbsStart);
    parm[7]  = rAbs;
    parm[8]  = zNose * TMath::Tan(theta1);

    parm[9]  = parm[3]+(zConeTPC-zAbsStart);
    parm[10] = rAbs;
    parm[11] = parm[8] + (parm[9] - parm[6]) * TMath::Tan(theta2);

    parm[12] = parm[3]+(zOpen-zAbsStart);
    parm[13] = rAbs;
    parm[14] = parm[11] + (parm[12] - parm[9]) * TMath::Tan(accMax);

    parm[15] = parm[3]+(zRear-dRear-zAbsStart);
    parm[16] = rAbs   + (parm[15] - parm[12]) * TMath::Tan(thetaOpen1) ;
    parm[17] = parm[14] + (parm[15] - parm[12]) * TMath::Tan(accMax);

    parm[18] = parm[3]+(zRear-dRear-zAbsStart);
    parm[19] = (zRear-dRear) * TMath::Tan(accMin);
    parm[20] = parm[14] + (parm[18] - parm[12]) * TMath::Tan(accMax);

    parm[21] = -parm[3];
    parm[22] = zRear* TMath::Tan(accMin);
    parm[23] = parm[20] + (parm[21] - parm[18]) * TMath::Tan(accMax);
  
    printf(" zRear, zAbsStart,  rAbs =  %f , %f , %f \n\n", 
             zRear, zAbsStart,  rAbs ); 
	      
   
    TPCON *abs0 = new TPCON("abs0", "abs0", "void", parm[0], parm[1], parm[2]);
    
    
    abs0->DefineSection(0, parm[3],  parm[4],  parm[5] );
    abs0->DefineSection(1, parm[6],  parm[7],  parm[8] );
    abs0->DefineSection(2, parm[9], parm[10], parm[11] );
    abs0->DefineSection(3,parm[12], parm[13], parm[14] );
    abs0->DefineSection(4,parm[15], parm[16], parm[17] );
    abs0->DefineSection(5,parm[18], parm[19], parm[20] );
    abs0->DefineSection(6,parm[21], parm[22], parm[23] );

    dz = (zRear-zAbsStart)/2.+zAbsStart;
    
    TRotMatrix* mat921 = new TRotMatrix("rot921","rot921",90.0,0.0,90.,90.0,180.0,0.0);
    
    Node0 = new TNode("abs0","abs0",abs0,0.0,0.0,dz,0);
    Node0->SetLineColor(38);	
    fNodes->Add(Node0);
    
    
    Float_t cpar[5];
    
    cpar[0] = (zNose - zAbsStart) / 2.;
    cpar[1] = zAbsStart * TMath::Tan(accMax);
    cpar[2] = zAbsStart * TMath::Tan(theta1)-dSteel;
    cpar[3] = zNose * TMath::Tan(accMax);
    cpar[4] = zNose * TMath::Tan(theta1)-dSteel;    

    dz = -(zRear-zAbsStart)/2.+cpar[0];
    
    TCONE *abs1 =  new TCONE("abs1", "abs1", "void", cpar[0], cpar[1], cpar[2], 
                              cpar[3], cpar[4]);        
    
    Node0->cd();
  
    Node1 = new TNode("abs1","abs1",abs1,0.0,0.0,dz,0);
    Node1->SetLineColor(7);	
    fNodes->Add(Node1);
    
//  Here I add a reference box to visualise the vertex zone :
    
    Top->cd();
    
    Float_t  paref[3];
    
    paref[0]  =  10.0;
    paref[1]  =  15.0; 
    paref[2]  =   5.6;
    
    TBRIK *aref = new TBRIK("aref","aref", "void", paref[0],paref[1],paref[2]);
			    
    Node2 = new TNode("aref","aref",aref,0.0,0.0,0.0,0);
    Node2->SetLineColor(kBlue);	
    fNodes->Add(Node2);    
      
//  Here I add the  mother volume  QBPM
//  and the flanges  QB29, QB22 et QB24 defined by Andreas Morsch 
//  in   AliPIPEv0::CreateGeometry()  :
 

    Float_t parp[36];
    
//  Mother Volume QBPM :

    parp[0]  =    0;
    parp[1]  =  360;
    parp[2]  =   11;
 
    parp[3]  = - 90;
    parp[4]  =    0;
    parp[5]  =    5.8;

    parp[6]  = - 81.0;
    parp[7]  =    0.;
    parp[8]  =    5.8;

    parp[9]  = - 81.;
    parp[10]  =    0.;
    parp[11] =    4.22;

    parp[12] = - 28.00;
    parp[13] =    0;
    parp[14] =    4.22;

    parp[15] = - 28.00;
    parp[16] =    0;
    parp[17] =    3.2;

    parp[18] =    0;
    parp[19] =    0;
    parp[20] =    3.2;

    parp[21] =   28.00;
    parp[22] =    0;
    parp[23] =    3.2;

    parp[24] =   28.00;
    parp[25] =    0;
    parp[26] =    4.22;

    parp[27] =  250;
    parp[28] =    0;
    parp[29] =   4.22;

    parp[30] =  250;
    parp[31] =    0;
    parp[32] =    5;

    parp[33] =  800;
    parp[34] =    0;
    parp[35] =    5;
    
    TPCON *pip0 = new TPCON("pip0", "pip0", "void", parp[0], parp[1], parp[2]);
     
    pip0->DefineSection( 0, parp[3],  parp[4],  parp[5] );
    pip0->DefineSection( 1, parp[6],  parp[7],  parp[8] );
    pip0->DefineSection( 2, parp[9], parp[10], parp[11] );
    pip0->DefineSection( 3,parp[12], parp[13], parp[14] );
    pip0->DefineSection( 4,parp[15], parp[16], parp[17] );
    pip0->DefineSection( 5,parp[18], parp[19], parp[20] );
    pip0->DefineSection( 6,parp[21], parp[22], parp[23] );
    pip0->DefineSection( 7,parp[24], parp[25], parp[26] );
    pip0->DefineSection( 8,parp[27], parp[28], parp[29] );
    pip0->DefineSection( 9,parp[30], parp[31], parp[32] );
    pip0->DefineSection(10,parp[33], parp[34], parp[35] );
    
    dz = 0.0;
    
    Top->cd();
    
    Node3 = new TNode("pip0","pip0",pip0,0.0,0.0,dz,mat921);
    Node3->SetLineColor(10);	
    fNodes->Add(Node3);
    Node3->SetVisibility(2);
    
//  Flanges QB29 at  654.8  and   254.8  cms :

    Float_t  ptube[3];
    
    ptube[0] = 3.0;
    ptube[1] = 4.9;
    ptube[2] = 2.2; 
    
    TTUBE *pip1 = new TTUBE("pip1", "pip1", "void", ptube[0], ptube[1], ptube[2]);  
     
    Node3->cd();
       
    
    Node4 = new TNode("pip1","pip1",pip1,0.0,0.0,254.8,0);
    Node4->SetLineColor(6);	
    fNodes->Add(Node4);
    
    TTUBE *pip2 = new TTUBE("pip2", "pip2", "void", ptube[0], ptube[1], ptube[2]);  
    
    Node3->cd();
    
    Node5 = new TNode("pip2","pip2",pip2,0.0,0.0,654.8,0);
    Node5->SetLineColor(6);	
    fNodes->Add(Node5);
    

//  Al-Be  section QBAB at 335.0 cm  (LEFT side) :
    
    ptube[0] =  2.90;
    ptube[1] =  3.05;
    ptube[2] = 171.0;

    TTUBE *pip3 = new TTUBE("pip3", "pip3", "void", ptube[0], ptube[1], ptube[2]); 
     
    Node3->cd();

    Node6 = new TNode("pip3","pip3",pip3,0.0,0.0,335.0+ptube[2],0);
    Node6->SetLineColor(6);
    fNodes->Add(Node6);
    
// Here I add the flange which is sitting on beam line 
// right in front of V0R detector, and which I found on CERN drawing 
// entitled : ALICE BEAM VACCUM CHAMBER - RB26 version III :  

    ptube[0] = 3.0;
    ptube[1] = 5.675;
    ptube[2] = 0.9; 
    
    TTUBE *pip4 = new TTUBE("pip4", "pip4", "void", ptube[0], ptube[1], ptube[2]);  
     
    Node3->cd();
    
    
    Node7 = new TNode("pip4","pip4",pip4,0.0,0.0,-85.0-0.9,0);
    Node7->SetLineColor(6);
    fNodes->Add(Node7);   
         
 }  
    

//------------------------------------------------------------------------
void AliVZEROv0::CreateMaterials()
{
  Int_t i;

  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" VZERO create materials ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");

/* ...................  OLD VALUES  ( used in RICH) ............................   
    Float_t ppckov[14] = { 5.63e-9,5.77e-9,5.9e-9,6.05e-9,6.2e-9,6.36e-9,6.52e-9,
			   6.7e-9,6.88e-9,7.08e-9,7.3e-9,7.51e-9,7.74e-9,8e-9 };			   
    Float_t rindex_quarz[14] = { 1.528309,1.533333,
				 1.538243,1.544223,1.550568,1.55777,
				 1.565463,1.574765,1.584831,1.597027,
			         1.611858,1.6277,1.6472,1.6724 };
    Float_t absco_quarz[14] = { 20.126,16.27,13.49,11.728,9.224,8.38,7.44,7.17,
				6.324,4.483,1.6,.323,.073,0. };	
...................................................................................... */				

             
    Float_t ppckov[14] = { 5.5e-9, 5.7e-9, 5.9e-9, 6.1e-9, 6.3e-9, 6.5e-9, 6.7e-9, 
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
    
    TGeant3 *geant3 = (TGeant3*) gMC;
    
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
       
    AliMaterial( 1, "AIR A$", 14.61, 7.3, .001205, 30420., 67500);
    AliMaterial(11, "AIR I$", 14.61, 7.3, .001205, 30420., 67500);
    AliMaterial( 2, "CARBON$"  , 12.01, 6.0, 2.265, 18.8, 49.9);
    AliMixture(  3, "QUA", aqua, zqua, densqua, nlmatqua, wmatqua);
    AliMaterial( 4, "ALUMINIUM1$", 26.98, 13., 2.7, 8.9, 37.2);
    AliMaterial( 5, "ALUMINIUM2$", aal, zal, densal, radlal, 0);
    
    
    AliMixture( 6, "Scintillator$",ascin,zscin,denscin,-2,wscin);
    
     
    Int_t   ISXFLD = gAlice->Field()->Integ();
    Float_t SXMGMX = gAlice->Field()->Max();
    
    Float_t tmaxfd, stemax, deemax, epsil, stmin;
    
    tmaxfd = -10.;
    stemax = -0.1;
    deemax = -0.1;
    epsil  = -0.01;
    stmin  = -0.001;
    

//  Active Air :    
    AliMedium(1, "ACTIVE AIR$", 1, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);

//  Inactive air : 
  
    AliMedium(11, "INACTIVE AIR$", 11, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
    
    AliMedium(2, "CARBON$ ", 2,  1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);   

    AliMedium(3, "QUARZ$", 3, 1, ISXFLD, SXMGMX, tmaxfd, fMaxStepQua, fMaxDestepQua, epsil, stmin);
    
    AliMedium(4,"ALUMINUM1$",4, 1, ISXFLD, SXMGMX, tmaxfd, fMaxStepAlu, fMaxDestepAlu, epsil, stmin);
    AliMedium(5,"ALUMINUM2$",5, 1, ISXFLD, SXMGMX, tmaxfd, fMaxStepAlu, fMaxDestepAlu, epsil, stmin);    
    AliMedium(6,"SCINTILLATOR$",6, 1, ISXFLD, SXMGMX, 10., .01, 1., .003, .003);

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
    
    gMC->Gstpar(idtmed[3003], "LOSS", 1.);  //  [3003] = normal aluminum
    gMC->Gstpar(idtmed[3003], "HADR", 1.);
    gMC->Gstpar(idtmed[3003], "DCAY", 1.);
    gMC->Gstpar(idtmed[3003], "DRAY", 1.);
    
    gMC->Gstpar(idtmed[3004], "LOSS", 1.);  //  [3004] = reflecting aluminum
    gMC->Gstpar(idtmed[3004], "HADR", 1.);
    gMC->Gstpar(idtmed[3004], "DCAY", 1.);
    gMC->Gstpar(idtmed[3004], "DRAY", 1.);
    
    gMC->Gstpar(idtmed[3005], "LOSS", 1.);  //  [3005] = scintillator
    gMC->Gstpar(idtmed[3005], "HADR", 1.);
    gMC->Gstpar(idtmed[3005], "DCAY", 1.);
    gMC->Gstpar(idtmed[3005], "DRAY", 1.);    
    
    geant3->Gsckov(idtmed[3002], 14, ppckov, absco_quarz, effic_all,rindex_quarz);
    
    geant3->Gsckov(idtmed[3004], 14, ppckov, absco_alu, effic_alu, rindex_alu);

}
//---------------------------------------------------------------------
void AliVZEROv0::DrawModule()
{

//  Drawing is done in DrawVZERO.C

   Int_t i;

   printf("\n");
   for(i=0;i<35;i++) printf("*");
   printf(" VZERO DrawModule ");
   for(i=0;i<35;i++) printf("*");
   printf("\n");


}

//-------------------------------------------------------------------
void AliVZEROv0::Init()
{
// Initialises version 0 of the VZERO Detector
// Just prints an information message
  
   Int_t i;

   printf("\n");
   for(i=0;i<35;i++) printf("*");
   printf(" VZERO_Init \n");
   for(i=0;i<35;i++) printf("*");
   printf("\n");
  
   fMulti       = 0;
   fNCerenkovs  = 0;
   fNGCerenkovs = 0;
   fNdead       = 0;
  
   BookingHistograms();
   
}

//-------------------------------------------------------------------

void AliVZEROv0::StepManager()
{
  
//   Minimal version of StepManager :
//   Everything has been removed, I only AddHit whenever hit is in 
//   volume V0RI.  
  
     Int_t          copy;
     Int_t          vol[4];  // (box, layer, row, column) indices
     Float_t        hits[19];  // position wrt MRS,   energies...
     
     TLorentzVector pos;
     Float_t        global[3];
     Float_t        local[3];
     
     TLorentzVector momentum;
     Float_t        theta;
     Float_t        phi;
     Float_t        mom[4];
     Float_t        kRaddeg = 180/TMath::Pi();
     
     Float_t        TrackEnters = 0.0;
     Float_t        TrackExits  = 0.0;
     Float_t        Cerenkov    = 0.0;
     
     gMC->SetMaxStep(fMaxStepAlu);
     gMC->SetMaxStep(fMaxStepQua);
     
     if (!gMC->IsTrackAlive()) return;
     
     if (gMC->IsTrackEntering())       TrackEnters = 1.0;
     if (gMC->IsTrackExiting() )       TrackExits  = 1.0;
     if (gMC->TrackPid() == 50000050)  Cerenkov    = 1.0;
     
  
        
     gMC->TrackPosition(pos);
     gMC->TrackMomentum(momentum);
     
     mom[0] = momentum[0];
     mom[1] = momentum[1];
     mom[2] = momentum[2];
     mom[3] = momentum[3];
     
     Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
     Double_t rt = TMath::Sqrt(tc);
     theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
     phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
    
     global[0] = pos[0];
     global[1] = pos[1];
     global[2] = pos[2];
     
     
     gMC->Gmtod(global,local,1);
     
     hits[0]   = global[0];
     hits[1]   = global[1];
     hits[2]   = global[2];
     hits[3]   = local[0];
     hits[4]   = local[1];
     hits[5]   = local[2];
     hits[6]   = gMC->Edep();
     hits[7]   = gMC->Etot();
     hits[8]   = Float_t (gMC->TrackPid());
     hits[9]   = Float_t (gMC->IdFromPDG(gMC->TrackPid()));
     hits[10]  = gMC->TrackTime();
     hits[11]  = TrackEnters;
     hits[12]  = TrackExits;
     hits[13]  = gMC->TrackCharge();
     hits[14]  = Cerenkov;  
     
     hits[16]  = theta;
     hits[17]  = phi; 
   
     vol[0]    = gMC->CurrentVolOffID(1, vol[1]);
     vol[2]    = gMC->CurrentVolID(copy);
     vol[3]    = copy;


     if (gMC->CurrentVolID(copy) >= gMC->VolId("V0RI") &&
         gMC->CurrentVolID(copy) <= gMC->VolId("V0E4"))	 
     {
    	 AddHit(gAlice->CurrentTrack(), vol, hits);
     	}       

  }
 
//---------------------------------------------------------------------
void AliVZEROv0::AddHit(Int_t track, Int_t* vol, Float_t* hits)
{
  
  // Adds a  hit 
  
  
   TClonesArray  &lhits = *fHits;
   
   new(lhits[fNhits++]) AliVZEROhit(fIshunt, track, vol, hits);

}

//---------------------------------------------------------------------
void AliVZEROv0::FinishEvent()
{
   
   printf("\n");
   for(int i=0;i<30;i++) printf("*");
   printf(" VZERO_finishevent");
   for(int i=0;i<30;i++) printf("*");
   printf("\n");
   
     AddDigit(tracks, digits);
     
     
     if(fMulti > 0) fhMultiplicity->Fill(fMulti); 
     fhCerenkov->Fill(fNCerenkovs); 

     fMulti       = 0;
     fNCerenkovs  = 0;
     fNGCerenkovs = 0;
     fNdead       = 0;
     
}

//---------------------------------------------------------------------
void AliVZEROv0::AddDigit(Int_t *tracks, Int_t* digits) 
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
  
  
  char *H = strstr(option,"H");
  
  if (fHits   && gAlice->TreeH() && H) {
    gAlice->TreeH()->Branch(branchname,&fHits, fBufferSize);
    printf("* AliDetector::MakeBranch * Making Branch %s for hits\n",branchname);
  }     

  char *D = strstr(option,"D");
  //
  if (fDigits   && gAlice->TreeD() && D) {
    gAlice->TreeD()->Branch(branchname,&fDigits, fBufferSize);
    printf("* AliDetector::MakeBranch * Making Branch %s for digits\n",branchname);
  }     


}

//---------------------------------------------------------------------
void AliVZEROv0::BookingHistograms()
{

  if (fhMultiplicity)         
    delete fhMultiplicity;
    
  if (fhGEANTcode) 
    delete fhGEANTcode;
    
  if (fhCerenkov)
    delete fhCerenkov; 
     
  if (fhToF)
    delete fhToF;    

//  fhMultiplicity = new TH1F("hMultiplicity", "hMultiplicity", 350,  0. , 350.); 

  fhMultiplicity = new TH1F("hMultiplicity", "hMultiplicity", 100,  1. , 100.);
  fhGEANTcode    = new TH1F("hGEANTcode", "hGEANTcode", 50,  1., 50.);
  fhCerenkov     = new TH1F("hCerenkov", "hCerenkov", 100, 1., 100000.);
  fhToF          = new TH1F("hToF", "hToF",150,2.0,7.0);
  
}
  
//---------------------------------------------------------------------
void AliVZEROv0::FinishRun()
{

   SavingHistograms();
}



//---------------------------------------------------------------------
void AliVZEROv0::SavingHistograms()
{

//  Saves the histograms in a root file named "name.save" 


  Text_t outputname[8] ;
  outputname = "Fileout";
  TFile output(outputname,"RECREATE");
  output.cd();
  
  if (fhMultiplicity)         
    fhMultiplicity->Write();
  if (fhGEANTcode)
    fhGEANTcode->Write(); 
  if (fhCerenkov)
    fhCerenkov->Write();
  if (fhToF)
    fhToF->Write();   
}
