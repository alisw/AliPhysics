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
/*
$Log$
Revision 1.27  2003/10/08 12:56:58  bnandi
gaspmd[2] value changed from 7.0 to 0.25

Revision 1.26  2003/10/03 06:04:10  bnandi
z_psa and z_psb bugs fixed

Revision 1.25  2003/10/01 11:08:04  bnandi
changes for NewIO

Revision 1.24  2003/10/01 08:32:51  hristov
CurrentTrack replaced by GetCurrentTrackNumber

Revision 1.23  2003/10/01 05:07:51  bnandi
New geometry in new Alice Coordinate system

New rectangular geometry for ALICE PMD - Bedanga Mohanty and Y. P. Viyogi
June 2003
*/
//
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Photon Multiplicity Detector Version 1                                   //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliPMDv1Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
////

#include "AliPMDv1.h"
#include "AliRun.h"
#include "AliConst.h" 
#include "AliMagF.h" 
#include "Riostream.h"
#include <TVirtualMC.h>
 
static Int_t     ncol_um1,ncol_um2, nrow_um1, nrow_um2;
static Int_t     kdet;
static Float_t   sm_length_ax,sm_length_ay;
static Float_t   sm_length_bx,sm_length_by;
static Float_t   zdist, zdist1;
static Float_t   sm_thick, cell_radius, cell_wall, cell_depth;
static Float_t   boundary, th_base, th_air, th_pcb;
static Float_t   th_lead, th_steel;

ClassImp(AliPMDv1)
 
  //_____________________________________________________________________________
  AliPMDv1::AliPMDv1()
{
  //
  // Default constructor 
  //
  fMedSens=0;
}
 
//_____________________________________________________________________________
AliPMDv1::AliPMDv1(const char *name, const char *title)
  : AliPMD(name,title)
{
  //
  // Standard constructor
  //
  fMedSens=0;
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
  //components. They have 9 unit moudles inside them
  // ESMA, ESMB are placed in EPMD along with EMPB (Pb converter) 
  // and EMFE (iron support) 

  
  Int_t i,j;
  Float_t xb, yb, zb;
  Int_t number;
  Int_t ihrotm,irotdm;
  const Float_t root3_2 = TMath::Sqrt(3.) /2.; 
  const Float_t root3 = TMath::Sqrt(3.); 
  Int_t *idtmed = fIdtmed->GetArray()-599;
 
  AliMatrix(ihrotm, 90., 30.,   90.,  120., 0., 0.);
  AliMatrix(irotdm, 90., 180.,  90.,  270., 180., 0.);
 
  zdist = TMath::Abs(zdist1);

  // First create the sensitive medium of a hexagon cell (ECAR)
  // Inner hexagon filled with gas (Ar+CO2)
  
  Float_t hexd2[10] = {0.,360.,6,2,-0.25,0.,0.23,0.25,0.,0.23};
  hexd2[4] = -cell_depth/2.;
  hexd2[7] =  cell_depth/2.;
  hexd2[6] =  cell_radius - cell_wall;
  hexd2[9] =  cell_radius - cell_wall;
  
  gMC->Gsvolu("ECAR", "PGON", idtmed[604], hexd2,10);
  gMC->Gsatt("ECAR", "SEEN", 0);
  
  // Place the sensitive medium inside a hexagon copper cell (ECCU)
  // Outer hexagon made of Copper
  
  Float_t hexd1[10] = {0.,360.,6,2,-0.25,0.,0.25,0.25,0.,0.25};
  hexd1[4] = -cell_depth/2.;
  hexd1[7] =  cell_depth/2.;
  hexd1[6] =  cell_radius;
  hexd1[9] =  cell_radius;

  gMC->Gsvolu("ECCU", "PGON", idtmed[614], hexd1,10);
  gMC->Gsatt("ECCU", "SEEN", 0);

  // Place  inner hex (sensitive volume) inside outer hex (copper)
  
  gMC->Gsposp("ECAR", 1, "ECCU", 0., 0., 0., 0, "ONLY", hexd2, 10);
  
  // Now create Rectangular TWO strips (EST1, EST2) 
  // of 1 column and 48 or 96 cells length

  // volume for first strip EST1 made of AIR 

  Float_t dbox1[3];
  dbox1[0] = ncol_um1*cell_radius;
  dbox1[1] = cell_radius/root3_2;
  dbox1[2] = cell_depth/2.;
  
  gMC->Gsvolu("EST1","BOX", idtmed[698], dbox1, 3);
  gMC->Gsatt("EST1", "SEEN", 0);

  // volume for second strip EST2 

  Float_t dbox2[3];
  dbox2[0] = ncol_um2*cell_radius;
  dbox2[1] = dbox1[1];
  dbox2[2] = dbox1[2];

  gMC->Gsvolu("EST2","BOX", idtmed[698], dbox2, 3);
  gMC->Gsatt("EST2", "SEEN", 0);

  // Place hexagonal cells ECCU placed inside EST1 
  yb = 0.; 
  zb = 0.;
  xb = -(dbox1[0]) + cell_radius; 
  for (i = 1; i <= ncol_um1; ++i) 
	{
	  number = i;
	  gMC->Gsposp("ECCU", number, "EST1", xb,yb,zb, ihrotm, "ONLY", hexd1,10);
	  xb += (cell_radius*2.);
	}
  // Place hexagonal cells ECCU placed inside EST2 
      yb = 0.; 
      zb = 0.;
      xb = -(dbox2[0]) + cell_radius; 
      for (i = 1; i <= ncol_um2; ++i) 
	{
	  number = i;
	  gMC->Gsposp("ECCU", number, "EST2", xb,yb,zb, ihrotm, "ONLY", hexd1,10);
	  xb += (cell_radius*2.);
	}



  // 2 types of rectangular shaped unit modules EUM1 and EUM2 (defined by BOX) 

  // Create EUM1

  Float_t dbox3[3];
  dbox3[0] = dbox1[0]+cell_radius/2.;
  dbox3[1] = (dbox1[1]*nrow_um1)-(cell_radius*root3*(nrow_um1-1)/6.); 
  dbox3[2] = cell_depth/2.;
  
  gMC->Gsvolu("EUM1","BOX", idtmed[698], dbox3, 3);
  gMC->Gsatt("EUM1", "SEEN", 1);
  
  // Place rectangular strips EST1 inside EUM1 unit module

  yb = -dbox3[1]+dbox1[1];  
  for (j = 1; j <= nrow_um1; ++j)  
    {
      if(j%2 == 0)
	{
      xb =cell_radius/2.0;
	}
      else
	{
	  xb = -cell_radius/2.0;
	}
      number = j;
      gMC->Gsposp("EST1",number, "EUM1", xb, yb , 0. , 0, "MANY",dbox1,3);
      yb = (-dbox3[1]+dbox1[1])+j*1.0*cell_radius*root3;
    }

  // Create EUM2

  Float_t dbox4[3];
  dbox4[0] = dbox2[0]+cell_radius/2.;
  dbox4[1] =(dbox2[1]*nrow_um2)-(cell_radius*root3*(nrow_um2-1)/6.); 
  dbox4[2] = dbox3[2];

  gMC->Gsvolu("EUM2","BOX", idtmed[698], dbox4, 3);
  gMC->Gsatt("EUM2", "SEEN", 1);

  // Place rectangular strips EST2 inside EUM2 unit module

  yb = -dbox4[1]+dbox2[1]; 
  for (j = 1; j <= nrow_um2; ++j) 
      {
      if(j%2 == 0)
	{
      xb =cell_radius/2.0;
	}
      else
	{
	  xb = -cell_radius/2.0;
	}
      number = j;
      gMC->Gsposp("EST2",number, "EUM2", xb, yb , 0. , 0, "MANY",dbox2,3);
      yb = (-dbox4[1]+dbox2[1])+j*1.0*cell_radius*root3;
    }

  // 2 types of Rectangular shaped supermodules (BOX) 
  //each with 6 unit modules 
  
  // volume for SUPERMODULE ESMA 
  //Space added to provide a gapping for HV between UM's

  Float_t dbox_sm1[3];
  dbox_sm1[0] = 3.0*dbox3[0]+(2.0*0.025);
  dbox_sm1[1] = 2.0*dbox3[1]+0.025;
  dbox_sm1[2] = cell_depth/2.;

  gMC->Gsvolu("ESMA","BOX", idtmed[698], dbox_sm1, 3);
  gMC->Gsatt("ESMA", "SEEN", 1);

  //Position the 6 unit modules in EMSA
  Float_t x_a1,x_a2,x_a3,y_a1,y_a2; 
  x_a1 = -dbox_sm1[0] + dbox3[0];
  x_a2 = 0.;
  x_a3 = dbox_sm1[0]  - dbox3[0]; 
  y_a1 = dbox_sm1[1]  - dbox3[1];
  y_a2 = -dbox_sm1[1] + dbox3[1];
  
  gMC->Gsposp("EUM1", 1, "ESMA", x_a1, y_a1, 0., 0, "ONLY",dbox3,3);
  gMC->Gsposp("EUM1", 2, "ESMA", x_a2, y_a1, 0., 0, "ONLY",dbox3,3);
  gMC->Gsposp("EUM1", 3, "ESMA", x_a3, y_a1, 0., 0, "ONLY",dbox3,3);
  gMC->Gsposp("EUM1", 4, "ESMA", x_a1, y_a2, 0., 0, "ONLY",dbox3,3);
  gMC->Gsposp("EUM1", 5, "ESMA", x_a2, y_a2, 0., 0, "ONLY",dbox3,3);
  gMC->Gsposp("EUM1", 6, "ESMA", x_a3, y_a2, 0., 0, "ONLY",dbox3,3);


  // volume for SUPERMODULE ESMB 
  //Space is added to provide a gapping for HV between UM's
  Float_t dbox_sm2[3];
  dbox_sm2[0] = 2.0*dbox4[0]+0.025;
  dbox_sm2[1] = 3.0*dbox4[1]+(2.0*0.025);
  dbox_sm2[2] = cell_depth/2.;
  
  gMC->Gsvolu("ESMB","BOX", idtmed[698], dbox_sm2, 3);
  gMC->Gsatt("ESMB", "SEEN", 1);

  //Position the 6 unit modules in EMSB
  Float_t x_b1,x_b2,y_b1,y_b2,y_b3; 
  x_b1 = -dbox_sm2[0] +dbox4[0];
  x_b2 = dbox_sm2[0]-dbox4[0];
  y_b1  =dbox_sm2[1]-dbox4[1];
  y_b2  = 0.; 
  y_b3  = -dbox_sm2[1]+dbox4[1];
  
  gMC->Gsposp("EUM2", 1, "ESMB", x_b1, y_b1, 0., 0, "ONLY",dbox4,3);
  gMC->Gsposp("EUM2", 2, "ESMB", x_b2, y_b1, 0., 0, "ONLY",dbox4,3);
  gMC->Gsposp("EUM2", 3, "ESMB", x_b1, y_b2, 0., 0, "ONLY",dbox4,3);
  gMC->Gsposp("EUM2", 4, "ESMB", x_b2, y_b2, 0., 0, "ONLY",dbox4,3);
  gMC->Gsposp("EUM2", 5, "ESMB", x_b1, y_b3, 0., 0, "ONLY",dbox4,3);
  gMC->Gsposp("EUM2", 6, "ESMB", x_b2, y_b3, 0., 0, "ONLY",dbox4,3);


  // Make a 3mm thick G10 Base plate for ESMA
  Float_t dbox_g1a[3];
  dbox_g1a[0] = dbox_sm1[0]; 
  dbox_g1a[1] = dbox_sm1[1];       
  dbox_g1a[2] = th_base/2.;

  gMC->Gsvolu("EBPA","BOX", idtmed[607], dbox_g1a, 3);
  gMC->Gsatt("EBPA", "SEEN", 1);

  // Make a 1.6mm thick G10 PCB for ESMA
  Float_t dbox_g2a[3];
  dbox_g2a[0] = dbox_sm1[0]; 
  dbox_g2a[1] = dbox_sm1[1];       
  dbox_g2a[2] = th_pcb/2.;

  gMC->Gsvolu("EPCA","BOX", idtmed[607], dbox_g2a, 3);
  gMC->Gsatt("EPCA", "SEEN", 1);


  // Make a Full module EFPA of AIR to place EBPA, 
  // 1mm AIR, EPCA, ESMA,EPCA for PMD
  
  Float_t dbox_alla[3];
  dbox_alla[0] = dbox_sm1[0]; 
  dbox_alla[1] = dbox_sm1[1];       
  dbox_alla[2] = (th_base+0.1+th_pcb+dbox_sm1[2]+th_pcb)/2.;

  gMC->Gsvolu("EFPA","BOX", idtmed[698], dbox_alla, 3);
  gMC->Gsatt("EFPA", "SEEN", 1);


  // Make a Full module EFCA of AIR to place EBPA, 
  // 1mm AIR, EPCA, ESMA,EPC for CPV
  Float_t dbox_alla2[3];
  dbox_alla2[0] = dbox_sm1[0]; 
  dbox_alla2[1] = dbox_sm1[1];       
  dbox_alla2[2] = (th_base+0.1+th_pcb+dbox_sm1[2]+th_pcb)/2.;

  gMC->Gsvolu("EFCA","BOX", idtmed[698], dbox_alla2, 3);
  gMC->Gsatt("EFCA", "SEEN", 1);

  // Now place everything in EFPA for PMD

  Float_t z_bpa,z_pcba1,z_pcba2,z_sma; 
  z_pcba1 = - dbox_alla[2]+th_pcb/2.0;
  gMC->Gsposp("EPCA", 1, "EFPA", 0., 0., z_pcba1, 0, "ONLY",dbox_g2a,3);
  z_sma = z_pcba1+dbox_sm1[2];
  gMC->Gsposp("ESMA", 1, "EFPA", 0., 0., z_sma, 0, "ONLY",dbox_sm1,3);
  z_pcba2 = z_sma+th_pcb/2.0;
  gMC->Gsposp("EPCA", 2, "EFPA", 0., 0., z_pcba2, 0, "ONLY",dbox_g2a,3);
  z_bpa = z_pcba2+0.1+th_base/2.0; // 0.1 for 0.1 mm Air gap 
  gMC->Gsposp("EBPA", 1, "EFPA", 0., 0., z_bpa, 0, "ONLY",dbox_g1a,3);

  // Now place everything in EFCA for CPV

  Float_t z_bpa2,z_pcba12,z_pcba22,z_sma2; 
  z_bpa2 = - dbox_alla2[2]+th_base/2.0;
  gMC->Gsposp("EBPA", 1, "EFCA", 0., 0., z_bpa2, 0, "ONLY",dbox_g1a,3);
  z_pcba12 = z_bpa2+0.1+th_pcb/2.0;
  gMC->Gsposp("EPCA", 1, "EFCA", 0., 0., z_pcba12, 0, "ONLY",dbox_g2a,3);
  z_sma2 = z_pcba12+dbox_sm1[2];
  gMC->Gsposp("ESMA", 1, "EFCA", 0., 0., z_sma2, 0, "ONLY",dbox_sm1,3);
  z_pcba22 = z_sma2+th_pcb/2.0;
  gMC->Gsposp("EPCA", 2, "EFCA", 0., 0., z_pcba22, 0, "ONLY",dbox_g2a,3);



  // Make a 3mm thick G10 Base plate for ESMB
  Float_t dbox_g1b[3];
  dbox_g1b[0] = dbox_sm2[0]; 
  dbox_g1b[1] = dbox_sm2[1];       
  dbox_g1b[2] = th_base/2.;

  gMC->Gsvolu("EBPB","BOX", idtmed[607], dbox_g1b, 3);
  gMC->Gsatt("EBPB", "SEEN", 1);

  // Make a 1.6mm thick G10 PCB for ESMB
  Float_t dbox_g2b[3];
  dbox_g2b[0] = dbox_sm2[0]; 
  dbox_g2b[1] = dbox_sm2[1];       
  dbox_g2b[2] = th_pcb/2.;

  gMC->Gsvolu("EPCB","BOX", idtmed[607], dbox_g2b, 3);
  gMC->Gsatt("EPCB", "SEEN", 1);


  // Make a Full module EFPB of AIR to place EBPB, 
  //1mm AIR, EPCB, ESMB,EPCB for PMD
  Float_t dbox_allb[3];
  dbox_allb[0] = dbox_sm2[0]; 
  dbox_allb[1] = dbox_sm2[1];       
  dbox_allb[2] = (th_base+0.1+th_pcb+dbox_sm2[2]+th_pcb)/2.;

  gMC->Gsvolu("EFPB","BOX", idtmed[698], dbox_allb, 3);
  gMC->Gsatt("EFPB", "SEEN", 1);

  // Make a Full module EFCB of AIR to place EBPB, 
  //1mm AIR, EPCB, ESMB,EPCB for CPV
  Float_t dbox_allb2[3];
  dbox_allb2[0] = dbox_sm2[0]; 
  dbox_allb2[1] = dbox_sm2[1];       
  dbox_allb2[2] = (th_base+0.1+th_pcb+dbox_sm2[2]+th_pcb)/2.;

  gMC->Gsvolu("EFCB","BOX", idtmed[698], dbox_allb2, 3);
  gMC->Gsatt("EFCB", "SEEN", 1);


  // Now place everything in EFPB for PMD

  Float_t z_bpb,z_pcbb1,z_pcbb2,z_smb; 
  z_pcbb1 = - dbox_allb[2]+th_pcb/2.0;
  gMC->Gsposp("EPCB", 1, "EFPB", 0., 0., z_pcbb1, 0, "ONLY",dbox_g2b,3);
  z_smb = z_pcbb1+dbox_sm2[2];
  gMC->Gsposp("ESMB", 1, "EFPB", 0., 0., z_smb, 0, "ONLY",dbox_sm2,3);
  z_pcbb2 = z_smb+th_pcb/2.0;
  gMC->Gsposp("EPCB", 2, "EFPB", 0., 0., z_pcbb2, 0, "ONLY",dbox_g2b,3);
  z_bpb = z_pcbb2+0.1+th_base/2.0; // 0.1 for 0.1 mm Air gap 
  gMC->Gsposp("EBPB", 1, "EFPB", 0., 0., z_bpb, 0, "ONLY",dbox_g1b,3);


  // Now place everything in EFCB for CPV

  Float_t z_bpb2,z_pcbb12,z_pcbb22,z_smb2; 
  z_bpb2 = - dbox_allb2[2]+th_base/2.0;
  gMC->Gsposp("EBPB", 1, "EFCB", 0., 0., z_bpb2, 0, "ONLY",dbox_g1b,3);
  z_pcbb12 = z_bpb2+0.1+th_pcb/2.0;
  gMC->Gsposp("EPCB", 1, "EFCB", 0., 0., z_pcbb12, 0, "ONLY",dbox_g2b,3);
  z_smb2 = z_pcbb12+dbox_sm2[2];
  gMC->Gsposp("ESMB", 1, "EFCB", 0., 0., z_smb2, 0, "ONLY",dbox_sm2,3);
  z_pcbb22 = z_smb2+th_pcb/2.0;
  gMC->Gsposp("EPCB", 2, "EFCB", 0., 0., z_pcbb22, 0, "ONLY",dbox_g2b,3);


  // Master MODULE EMPA of aluminum for PMD
  //Float_t dbox_mm1[3];
  dbox_mm1[0] = dbox_sm1[0]+boundary; 
  dbox_mm1[1] = dbox_sm1[1]+boundary;       
  dbox_mm1[2] = dbox_alla[2];

  gMC->Gsvolu("EMPA","BOX", idtmed[603], dbox_mm1, 3);
  gMC->Gsatt("EMPA", "SEEN", 1);

  // Master MODULE EMCA of aluminum for CPV
  //Float_t dbox_mm12[3];
  dbox_mm12[0] = dbox_sm1[0]+boundary; 
  dbox_mm12[1] = dbox_sm1[1]+boundary;       
  dbox_mm12[2] = dbox_alla[2];

  gMC->Gsvolu("EMCA","BOX", idtmed[603], dbox_mm12, 3);
  gMC->Gsatt("EMCA", "SEEN", 1);


  //Position EFMA inside EMMA for PMD and CPV
  gMC->Gsposp("EFPA", 1, "EMPA", 0., 0., 0., 0, "ONLY",dbox_alla,3);
  gMC->Gsposp("EFCA", 1, "EMCA", 0., 0., 0., 0, "ONLY",dbox_alla2,3);


  // Master MODULE EMPB of aluminum for PMD
  //Float_t dbox_mm2[3];
  dbox_mm2[0] = dbox_sm2[0]+boundary; 
  dbox_mm2[1] = dbox_sm2[1]+boundary;       
  dbox_mm2[2] = dbox_allb[2];

  gMC->Gsvolu("EMPB","BOX", idtmed[603], dbox_mm2, 3);
  gMC->Gsatt("EMPB", "SEEN", 1);

  // Master MODULE EMCB of aluminum for CPV
  //Float_t dbox_mm22[3];
  dbox_mm22[0] = dbox_sm2[0]+boundary; 
  dbox_mm22[1] = dbox_sm2[1]+boundary;       
  dbox_mm22[2] = dbox_allb[2];

  gMC->Gsvolu("EMCB","BOX", idtmed[603], dbox_mm22, 3);
  gMC->Gsatt("EMCB", "SEEN", 1);

 
  //Position EFMB inside EMMB
  gMC->Gsposp("EFPB", 1, "EMPB", 0., 0., 0., 0, "ONLY",dbox_allb,3);
  gMC->Gsposp("EFCB", 1, "EMCB", 0., 0., 0., 0, "ONLY",dbox_allb2,3);

}
 
//_____________________________________________________________________________

void AliPMDv1::CreatePMD()
{

  //
  // Create final detector from supermodules
  // -- Author : Bedanga and Viyogi June 2003

  Float_t  xp, yp, zp;
  Int_t jhrot12,jhrot13, irotdm;
  Int_t *idtmed = fIdtmed->GetArray()-599;
  
  //VOLUMES Names : begining with "E" for all PMD volumes, 

  // --- DEFINE Iron, and lead volumes  for SM A
  
  Float_t dbox_pba[3];
  dbox_pba[0] = sm_length_ax;
  dbox_pba[1] = sm_length_ay;
  dbox_pba[2] = th_lead/2.;
  
  gMC->Gsvolu("EPBA","BOX", idtmed[600], dbox_pba, 3);
  gMC->Gsatt ("EPBA", "SEEN", 0);
  
  //   Fe Support 
  Float_t dbox_fea[3];
  dbox_fea[0] = sm_length_ax;
  dbox_fea[1] = sm_length_ay;
  dbox_fea[2] = th_steel/2.;
  
  gMC->Gsvolu("EFEA","BOX", idtmed[618], dbox_fea, 3);
  gMC->Gsatt ("EFEA", "SEEN", 0);

  // --- DEFINE Iron, and lead volumes  for SM B

  Float_t dbox_pbb[3];
  dbox_pbb[0] = sm_length_bx;
  dbox_pbb[1] = sm_length_by;
  dbox_pbb[2] = th_lead/2.;
  
  gMC->Gsvolu("EPBB","BOX", idtmed[600], dbox_pbb, 3);
  gMC->Gsatt ("EPBB", "SEEN", 0);
  
  //   Fe Support 
  Float_t dbox_feb[3];
  dbox_feb[0] = sm_length_bx;
  dbox_feb[1] = sm_length_by;
  dbox_feb[2] = th_steel/2.;
  
  gMC->Gsvolu("EFEB","BOX", idtmed[618], dbox_feb, 3);
  gMC->Gsatt ("EFEB", "SEEN", 0);


  // Gaspmd, the dimension of RECTANGULAR mother volume of PMD,

  //  Float_t gaspmd[3] = {81.5,94.5,7.};
  Float_t gaspmd[3] = {81.5,94.5,0.25};
  gaspmd[0] = sm_length_ax+sm_length_bx;
  gaspmd[1] = sm_length_ay+sm_length_by;


  gMC->Gsvolu("EPMD", "BOX", idtmed[698], gaspmd, 3);
  gMC->Gsatt("EPMD", "SEEN", 1);

  AliMatrix(irotdm, 90., 0.,  90.,  90., 180., 0.);
   
  AliMatrix(jhrot12, 90., 180., 90., 270., 0., 0.);
  AliMatrix(jhrot13, 90., 240., 90., 330., 0., 0.);

  Float_t x_sma,y_sma;
  Float_t x_smb,y_smb;
  x_sma = -(sm_length_bx)/1.0;
  y_sma = sm_length_by;
  x_smb = -sm_length_ax;
  y_smb = -sm_length_ay;

  //Complete detector for Type A
  //Position Super modules type A for both CPV and PMD in EPMD  
  Float_t z_psa,z_pba,z_fea,z_cva; 

  z_psa = - gaspmd[2] + sm_thick/2.;

  gMC->Gsposp("EMPA", 1, "EPMD", x_sma, y_sma, z_psa, 0, "ONLY",dbox_mm1,3);
  gMC->Gsposp("EMPA", 2, "EPMD", -x_sma, -y_sma, z_psa, jhrot12, "ONLY",dbox_mm1,3);
  z_pba=z_psa+sm_thick/2.+dbox_pba[2];
  gMC->Gsposp("EPBA", 1, "EPMD", x_sma, y_sma, z_pba, 0, "ONLY",dbox_pba,3);
  gMC->Gsposp("EPBA", 2, "EPMD", -x_sma, -y_sma, z_pba, 0, "ONLY",dbox_pba,3);
  z_fea=z_pba+dbox_pba[2]+dbox_fea[2];
  gMC->Gsposp("EFEA", 1, "EPMD", x_sma, y_sma, z_fea, 0, "ONLY",dbox_fea,3);
  gMC->Gsposp("EFEA", 2, "EPMD", -x_sma, -y_sma, z_fea, 0, "ONLY",dbox_fea,3);
  z_cva=z_fea+dbox_fea[2]+sm_thick/2.;
  gMC->Gsposp("EMCA", 1, "EPMD", x_sma, y_sma, z_cva, 0, "ONLY",dbox_mm12,3);
  gMC->Gsposp("EMCA", 2, "EPMD", -x_sma,-y_sma, z_cva, jhrot12, "ONLY",dbox_mm12,3);
 
  //Complete detector for Type B
  //Position Super modules type B for both CPV and PMD in EPMD  
  Float_t z_psb,z_pbb,z_feb,z_cvb; 
  z_psb = - gaspmd[2] + sm_thick/2.;
  
  gMC->Gsposp("EMPB", 3, "EPMD", x_smb, y_smb, z_psb, 0, "ONLY",dbox_mm2,3);
  gMC->Gsposp("EMPB", 4, "EPMD", -x_smb, -y_smb, z_psb, jhrot12, "ONLY",dbox_mm2,3);
  z_pbb=z_psb+sm_thick/2.+dbox_pbb[2];
  gMC->Gsposp("EPBB", 3, "EPMD", x_smb, y_smb, z_pbb, 0, "ONLY",dbox_pbb,3);
  gMC->Gsposp("EPBB", 4, "EPMD", -x_smb, -y_smb, z_pbb, 0, "ONLY",dbox_pbb,3);
  z_feb=z_pbb+dbox_pbb[2]+dbox_feb[2];
  gMC->Gsposp("EFEB", 3, "EPMD", x_smb, y_smb, z_feb, 0, "ONLY",dbox_feb,3);
  gMC->Gsposp("EFEB", 4, "EPMD", -x_smb, -y_smb, z_feb, 0, "ONLY",dbox_feb,3);
  z_cvb=z_feb+dbox_feb[2]+sm_thick/2.;
  gMC->Gsposp("EMCB", 3, "EPMD", x_smb, y_smb, z_cvb, 0, "ONLY",dbox_mm22,3);
  gMC->Gsposp("EMCB", 4, "EPMD", -x_smb,-y_smb, z_cvb, jhrot12, "ONLY",dbox_mm22,3);
  
  // --- Place the EPMD in ALICE 
  xp = 0.;
  yp = 0.;
  zp = zdist1;

  //Position Full PMD in ALICE   
  gMC->Gsposp("EPMD", 1, "ALIC", xp,yp,zp, 0, "ONLY",gaspmd,3);

}

 
//_____________________________________________________________________________
void AliPMDv1::DrawModule()
{
  cout << " Inside Draw Modules " << endl;
  //
  // Draw a shaded view of the Photon Multiplicity Detector
  //

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

  cout << " Outside Draw Modules " << endl;
}

//_____________________________________________________________________________
void AliPMDv1::CreateMaterials()
{
  cout << " Inside create materials " << endl;
  //
  // Create materials for the PMD
  //
  // ORIGIN    : Y. P. VIYOGI 
  //
  
  // --- The Argon- CO2 mixture --- 
  Float_t ag[2] = { 39.95 };
  Float_t zg[2] = { 18. };
  Float_t wg[2] = { .7,.3 };
  Float_t dar   = .001782;   // --- Ar density in g/cm3 --- 
  // --- CO2 --- 
  Float_t ac[2] = { 12.,16. };
  Float_t zc[2] = { 6.,8. };
  Float_t wc[2] = { 1.,2. };
  Float_t dc    = .001977;
  Float_t dco   = .002;  // --- CO2 density in g/cm3 ---
  
  Float_t absl, radl, a, d, z;
  Float_t dg;
  Float_t x0ar;
  Float_t buf[1];
  Int_t nbuf;
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  
  Int_t *idtmed = fIdtmed->GetArray()-599;
  Int_t isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  
  // --- Define the various materials for GEANT --- 
  AliMaterial(1, "Pb    $", 207.19, 82., 11.35, .56, 18.5);
  x0ar = 19.55 / dar;
  AliMaterial(2, "Argon$", 39.95, 18., dar, x0ar, 6.5e4);
  AliMixture(3, "CO2  $", ac, zc, dc, -2, wc);
  AliMaterial(4, "Al   $", 26.98, 13., 2.7, 8.9, 18.5);
  AliMaterial(6, "Fe   $", 55.85, 26., 7.87, 1.76, 18.5);
  AliMaterial(7, "W    $", 183.85, 74., 19.3, .35, 10.3);
  AliMaterial(8, "G10  $", 20., 10., 1.7, 19.4, 999.);
  AliMaterial(9, "SILIC$", 28.09, 14., 2.33, 9.36, 45.);
  AliMaterial(10, "Be   $", 9.01, 4., 1.848, 35.3, 36.7);
  AliMaterial(15, "Cu   $", 63.54, 29., 8.96, 1.43, 15.);
  AliMaterial(16, "C    $", 12.01, 6., 2.265, 18.8, 49.9);
  AliMaterial(17, "POLYCARBONATE    $", 20., 10., 1.2, 34.6, 999.);
  AliMixture(19, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel); 
  // AliMaterial(31, "Xenon$", 131.3, 54., dxe, x0xe, 6.5e4);
  
  AliMaterial(96, "MYLAR$", 8.73, 4.55, 1.39, 28.7, 62.);
  AliMaterial(97, "CONCR$", 20., 10., 2.5, 10.7, 40.);
  AliMaterial(98, "Vacum$", 1e-9, 1e-9, 1e-9, 1e16, 1e16);
  AliMaterial(99, "Air  $", 14.61, 7.3, .0012, 30420., 67500.);
 
  // 	define gas-mixtures 
  
  char namate[21];
  gMC->Gfmate((*fIdmate)[3], namate, a, z, d, radl, absl, buf, nbuf);
  ag[1] = a;
  zg[1] = z;
  dg = (dar * 4 + dco) / 5;
  AliMixture(5, "ArCO2$", ag, zg, dg, 2, wg);
  
  // Define tracking media 
  AliMedium(1, "Pb conv.$", 1,  0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
  AliMedium(7, "W  conv.$", 7,  0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
  AliMedium(8, "G10plate$", 8,  0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
  AliMedium(4, "Al      $", 4,  0, 0, isxfld, sxmgmx, .1,  .1, .01, .1);
  AliMedium(6, "Fe      $", 6,  0, 0, isxfld, sxmgmx, .1,  .1, .01, .1);
  AliMedium(5, "ArCO2   $", 5,  1, 0, isxfld, sxmgmx, .1,  .1, .1,  .1);
  AliMedium(9, "SILICON $", 9,  1, 0, isxfld, sxmgmx, .1,  .1, .1,  .1);
  AliMedium(10, "Be      $", 10, 0, 0, isxfld, sxmgmx, .1,  .1, .01, .1);
  AliMedium(98, "Vacuum  $", 98, 0, 0, isxfld, sxmgmx, 1., .1, .1,  10);
  AliMedium(99, "Air gaps$", 99, 0, 0, isxfld, sxmgmx, 1., .1, .1,  .1);
  AliMedium(15, "Cu      $", 15, 0, 0, isxfld, sxmgmx, .1,  .1, .01, .1);
  AliMedium(16, "C       $", 16, 0, 0, isxfld, sxmgmx, .1,  .1, .01, .1);
  AliMedium(17, "PLOYCARB$", 17, 0, 0, isxfld, sxmgmx, .1,  .1, .01, .1);
  AliMedium(19, " S steel$", 19, 0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
  //  AliMedium(31, "Xenon   $", 31,  1, 0, isxfld, sxmgmx, .1,  .1, .1,  .1);
  
  // --- Generate explicitly delta rays in the iron, aluminium and lead --- 
  gMC->Gstpar(idtmed[600], "LOSS", 3.);
  gMC->Gstpar(idtmed[600], "DRAY", 1.);
  
  gMC->Gstpar(idtmed[603], "LOSS", 3.);
  gMC->Gstpar(idtmed[603], "DRAY", 1.);
  
  gMC->Gstpar(idtmed[604], "LOSS", 3.);
  gMC->Gstpar(idtmed[604], "DRAY", 1.);
  
  gMC->Gstpar(idtmed[605], "LOSS", 3.);
  gMC->Gstpar(idtmed[605], "DRAY", 1.);
  
  gMC->Gstpar(idtmed[606], "LOSS", 3.);
  gMC->Gstpar(idtmed[606], "DRAY", 1.);
  
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
  gMC->Gstpar(idtmed[606], "CUTGAM", 1e-4);
  gMC->Gstpar(idtmed[606], "CUTELE", 1e-4);
  gMC->Gstpar(idtmed[606], "CUTNEU", 1e-4);
  gMC->Gstpar(idtmed[606], "CUTHAD", 1e-4);
  gMC->Gstpar(idtmed[603], "CUTGAM", 1e-4);
  gMC->Gstpar(idtmed[603], "CUTELE", 1e-4);
  gMC->Gstpar(idtmed[603], "CUTNEU", 1e-4);
  gMC->Gstpar(idtmed[603], "CUTHAD", 1e-4);
  gMC->Gstpar(idtmed[609], "CUTGAM", 1e-4);
  gMC->Gstpar(idtmed[609], "CUTELE", 1e-4);
  gMC->Gstpar(idtmed[609], "CUTNEU", 1e-4);
  gMC->Gstpar(idtmed[609], "CUTHAD", 1e-4);
  
  // --- Prevent particles stopping in the gas due to energy cut-off --- 
  gMC->Gstpar(idtmed[604], "CUTGAM", 1e-5);
  gMC->Gstpar(idtmed[604], "CUTELE", 1e-5);
  gMC->Gstpar(idtmed[604], "CUTNEU", 1e-5);
  gMC->Gstpar(idtmed[604], "CUTHAD", 1e-5);
  gMC->Gstpar(idtmed[604], "CUTMUO", 1e-5);

  cout << " Outside create materials " << endl;

}

//_____________________________________________________________________________
void AliPMDv1::Init()
{
  //
  // Initialises PMD detector after it has been built
  //

  Int_t i;
  kdet=1;
  //
  cout << " Inside Init " << endl;
  if(fDebug) {
      printf("\n%s: ",ClassName());
      for(i=0;i<35;i++) printf("*");
      printf(" PMD_INIT ");
      for(i=0;i<35;i++) printf("*");
      printf("\n%s: ",ClassName());
      printf("                 PMD simulation package (v1) initialised\n");
      printf("%s: parameters of pmd\n",ClassName());
      printf("%s: %10.2f %10.2f %10.2f \
      %10.2f\n",ClassName(),cell_radius,cell_wall,cell_depth,zdist1 );
      printf("%s: ",ClassName());
      for(i=0;i<80;i++) printf("*");
      printf("\n");
  }
  
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
  Int_t   vol[8]; //5
  //const char *namep;
  
  if(gMC->GetMedium() == fMedSens && (destep = gMC->Edep())) {
  
    gMC->CurrentVolID(copy);
    //namep=gMC->CurrentVolName();
    //printf("Current vol  is %s \n",namep);
    vol[0]=copy;

    gMC->CurrentVolOffID(1,copy);
    //namep=gMC->CurrentVolOffName(1);
    //printf("Current vol 11 is %s \n",namep);
    vol[1]=copy;

    gMC->CurrentVolOffID(2,copy);
    //namep=gMC->CurrentVolOffName(2);
    //printf("Current vol 22 is %s \n",namep);
    vol[2]=copy;

    //	if(strncmp(namep,"EHC1",4))vol[2]=1;

    gMC->CurrentVolOffID(3,copy);
    //namep=gMC->CurrentVolOffName(3);
    //printf("Current vol 33 is %s \n",namep);
    vol[3]=copy;

    gMC->CurrentVolOffID(4,copy);
    //namep=gMC->CurrentVolOffName(4);
    //printf("Current vol 44 is %s \n",namep);
    vol[4]=copy;

    gMC->CurrentVolOffID(5,copy);
    //namep=gMC->CurrentVolOffName(5);
    //printf("Current vol 55 is %s \n",namep);
    vol[5]=copy;

    gMC->CurrentVolOffID(6,copy);
    //namep=gMC->CurrentVolOffName(6);
    //printf("Current vol 66 is %s \n",namep);
    vol[6]=copy;

    gMC->CurrentVolOffID(7,copy);
    //namep=gMC->CurrentVolOffName(7);
    //printf("Current vol 77 is %s \n",namep);
    vol[7]=copy;


    //printf("volume number %4d %4d %4d %4d %4d %4d %4d %4d %10.3f \n",vol[0],vol[1],vol[2],vol[3],vol[4],vol[5],vol[6],vol[7],destep*1000000);
    
    gMC->Gdtom(center,hits,1);
    hits[3] = destep*1e9; //Number in eV
    AddHit(gAlice->GetCurrentTrackNumber(), vol, hits);

  }
}

  
//------------------------------------------------------------------------
// Get parameters

void AliPMDv1::GetParameters()
{
  const Float_t root3 = TMath::Sqrt(3.); 
  const Float_t root3_2 = TMath::Sqrt(3.) /2.; 
  //
  cell_radius=0.25;
  cell_wall=0.02;
  cell_depth=0.25 * 2.;
  //
  ncol_um1 = 48;
  ncol_um2 = 96;
  nrow_um1 = 96;//each strip has 1 row
  nrow_um2 = 48;//each strip has 1 row
  //
  sm_length_ax = (3.0*(ncol_um1*cell_radius+cell_radius/2.)+(2.0*0.025)) + 0.7;
  sm_length_bx = 2.0*(ncol_um2*cell_radius+cell_radius/2.)+0.025+0.7; 

  sm_length_ay = 2.0*(((cell_radius/root3_2)*nrow_um1)-(cell_radius*root3*(nrow_um1-1)/6.))+0.025+0.7;
  sm_length_by = 3.0*(((cell_radius/root3_2)*nrow_um2)-(cell_radius*root3*(nrow_um2-1)/6.))+(2.0*0.025)+0.7;
    //
    boundary=0.7;
    //
    th_base=0.3;
    th_air=0.1;
    th_pcb=0.16;
    //
    sm_thick = th_base + th_air + th_pcb + cell_depth + th_pcb + th_air + th_pcb;
    //
    th_lead=1.5;
    th_steel=0.5;
    
    zdist1 = 361.5;

}
