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

/*
$Log$
Revision 1.6  1999/09/29 09:24:28  fca
Introduction of the Copyright and cvs Log

*/

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

#include "AliPMDv1.h"
#include "AliRun.h"
#include "AliMC.h" 
#include "AliConst.h" 
 
static Int_t maxbox, kdet;
static Float_t thmin,thmax,zdist,zdist1,thlow,thhigh;

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
  //
  // Create geometry for Photon Multiplicity Detector Version 1
  //
  //Begin_Html
  /*
    <img src="picts/AliPMDv1.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliPMDv1Tree.gif">
  */
  //End_Html
  CreatePads();
  CreateInside();
}
 
//_____________________________________________________________________________
void AliPMDv1::CreateInside()
{
  //
  // Create inside of Pads
  //
  // -- Author :     Y.P. VIYOGI, 07/05/1996. 
  // -- Modified:    P.V.K.S.Baba(JU), 15-12-97. 
// Sipmd, the dimension of TUBE mother volume of PMD, other dimensions
// like sip01.. are to place more tubes in the volume at different eta bins.  
  Float_t sipmd[3] = { 40.,270.,15.};
  Float_t sip01[3] = { 10.,57.89,25.};
  Float_t sip02[3] = { 10.,64.03,25.};
  Float_t sip03[3] = { 10.,70.80,25.};
  Float_t sip04[3] = { 10.,78.32,25.};
  Float_t sip05[3] = { 10.,86.68,25.};
  Float_t sip06[3] = { 10.,95.91,25.};
  Float_t sip07[3] = { 10.,106.14,25.};
  Float_t sip08[3] = { 10.,117.48,25.};
  Float_t sip09[3] = { 10.,130.18,25.};
  Float_t sip10[3] = { 10.,144.18,25.};
  Float_t sip11[3] = { 10.,159.87,25.};
  Float_t sip12[3] = { 10.,177.43,25.};
  Float_t sip13[3] = { 10.,197.11,25.};
  Float_t sip14[3] = { 10.,219.28,25.};
  Float_t sipmdl[5] = { 10.,310.,25.,90.,270. };
  Float_t sipmdr[5] = { 10.,310.,25.,270.,90. };
  
  const Float_t root3_4 = sqrt(3)/4.;
  const Float_t root3_2 = sqrt(3)/2.;
  //  Float_t xiqa[4], yiqa[4];
  Int_t i;
  //  Float_t siqad[4];
  Float_t  xp, yp, zp;
  //  Int_t idrotm[100];
  Int_t num_mod;
  Int_t jhrotc,jhrotac;
//  const Float_t delx=78.8;
  const Float_t delx=76.75;
  //  const Float_t dely=delx*root3_2;
//  const Float_t delz=1.6/2.;
  AliMatrix(jhrotc, 90., 30.,   90.,  120., 0., 0.);
  AliMatrix(jhrotac, 90., 330., 90., 240., 0., 0.);
  Float_t x1= delx*root3_4;
  Float_t x2= delx*root3_4 + delx*root3_2;
  Float_t x3= delx*root3_4 + 2*delx*root3_2;
 Float_t xpos[13]={-x1,-x1,-x1,-x1,-x2,-x2,-x2,-x2,-x2,-x3,-x3,-x3,-x3};
  Float_t x4=delx/4.; 
 Float_t ypos[13]={(-70.-x4-delx),-(70.+x4),(70.+x4),(70.+x4+delx),-x4+2*delx,-x4+delx,-x4,-x4-delx,-x4-2*delx,-3*x4-delx,-x4-delx/2.,-3*x4+delx,-3*x4+2*delx};
// Float_t ypos[13]={(-70.-x4-delx),-(70.+x4),(70.+x4),(70.+x4+delx),(4*dely),(2*dely),0.,-(2*dely),-(4*dely),-3*x4-delx,-x4-delx/2.,-3*x4+delx,-3*x4+2*delx};
  Int_t *idtmed = fIdtmed->GetArray()-599;
  
  //  VOLUMES Names : begining with D for all PMD volumes, 
  // The names of SIZE variables begin with S and have more meaningful
  // characters as shown below. 
  
  // 		VOLUME 	SIZE	MEDIUM	: 	REMARKS 
  // 		------	-----	------	: --------------------------- 
  
  // 		DPMD	SIPMD	AIR	: INSIDE PMD  and its SIZE 
  
  
  
  // *** Define the  DPMD   Volume and fill with air *** 

  gMC->Gsvolu("DPMD", "TUBE", idtmed[698], sipmd, 3);
  gMC->Gsvolu("PM01", "TUBE", idtmed[698], sip01, 3);
  gMC->Gsvolu("PM02", "TUBE", idtmed[698], sip02, 3);
  gMC->Gsvolu("PM03", "TUBE", idtmed[698], sip03, 3);
  gMC->Gsvolu("PM04", "TUBE", idtmed[698], sip04, 3);
  gMC->Gsvolu("PM05", "TUBE", idtmed[698], sip05, 3);
  gMC->Gsvolu("PM06", "TUBE", idtmed[698], sip06, 3);
  gMC->Gsvolu("PM07", "TUBE", idtmed[698], sip07, 3);
  gMC->Gsvolu("PM08", "TUBE", idtmed[698], sip08, 3);
  gMC->Gsvolu("PM09", "TUBE", idtmed[698], sip09, 3);
  gMC->Gsvolu("PM10", "TUBE", idtmed[698], sip10, 3);
  gMC->Gsvolu("PM11", "TUBE", idtmed[698], sip11, 3);
  gMC->Gsvolu("PM12", "TUBE", idtmed[698], sip12, 3);
  gMC->Gsvolu("PM13", "TUBE", idtmed[698], sip13, 3);
  gMC->Gsvolu("PM14", "TUBE", idtmed[698], sip14, 3);
  gMC->Gsvolu("PMDL", "TUBS", idtmed[698], sipmdl, 5);
  gMC->Gsvolu("PMDR", "TUBS", idtmed[698], sipmdr, 5);
//  
  const Int_t npad2=72; 
  Float_t hexd1[10] = {0.,360.,6,2,-0.4,0.,0.53,0.4,0.,0.53};
  Float_t dpara_sm[6] = {12.5,12.5,0.8,30.,0.,0.};
  dpara_sm[0]=(npad2+0.25)*hexd1[6] + 1.2;
  dpara_sm[1] = dpara_sm[0] *root3_2;
  Float_t dpara_dm11[6] = {12.5,12.5,0.8,30.,0.,0.};
  dpara_dm11[0]=dpara_sm[0]+.01;
  dpara_dm11[1] = dpara_dm11[0] *root3_2;
  dpara_dm11[2]= 6.2/2.;
//
    for (i = 0; i < 2; ++i) {
        num_mod=i+1;
  gMC->Gsposp("DM11", num_mod, "DPMD", xpos[i],ypos[i],0., jhrotac, "ONLY", dpara_dm11, 6);
  gMC->Gsposp("DM11", num_mod+13, "DPMD", TMath::Abs(xpos[i]),ypos[i],0., jhrotc, "ONLY", dpara_dm11, 6);
    printf("Num_mod %d\n",num_mod);
	}
   maxbox=13;
    for (i = 2; i < maxbox; ++i) {
        num_mod=i+1;
  gMC->Gsposp("DM11", num_mod, "DPMD", xpos[i],ypos[i],0., jhrotc, "ONLY", dpara_dm11, 6);
  gMC->Gsposp("DM11", num_mod+13, "DPMD", TMath::Abs(xpos[i]),ypos[i],0., jhrotac, "ONLY", dpara_dm11, 6);
    printf("Num_mod %d\n",num_mod);
	}
//  gMC->Gspos("PM01", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM02", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM03", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM04", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM05", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM06", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM07", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM08", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM09", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM10", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM11", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM12", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM13", 1, "DPMD", 0.,0.,0., 0, "ONLY");
//  gMC->Gspos("PM14", 1, "DPMD", 0.,0.,0., 0, "ONLY");
// --- Place the DPMD in ALICE with front edge 5.8m from vertex  --- 
    xp = 0.;
    yp = 0.;
    zp = zdist1;
//  gMC->Gspos("PMDL", 1, "DPMD", xp,yp,0., 0, "ONLY");
//  gMC->Gspos("PMDR", 1, "DPMD", xp,yp,0., 0, "ONLY");
  gMC->Gspos("DPMD", 1, "ALIC", xp,yp,zp, 0, "ONLY");
    
}

//_____________________________________________________________________________
void AliPMDv1::CreatePads()
{
  //
  // Create the geometry of the pads
  // *** DEFINITION OF THE GEOMETRY OF THE PMD  *** 
  // *** HEXAGONAL PADS WITH 10 MM SQUARE EQUIVALENT
  // -- Author :     S. Chattopadhyay, 02/04/1999. 

// Basic unit is DP11, a hexagonal cell, which is placed inside another 
// hexagonal cell (DS11) of larger radius, compared to DP11. The difference in r// adius gives the dimension of half width of each cell wall.
// These cells are placed as 72 x 72 array in a 
// rhombus shaped supermodule (DW11). The rhombus shaped modules are designed
// to have closed packed structure.
// Each supermodule (SUPR), made of G10 is filled with following components
//  SMSS --> SS backing,
//  SMAR --> Gap between gas hexagonal cells and G10 backing.
//  DW11 --> Ar-Co2 filled gas hexagonal cells.
//  SMAR
// These supermodules are placed inside the main module (DM11), with Fe and 
// Pb converter positioned between CPV and PMD.
//  DM11 made of
// SUPR (rotated to place steel on the other side), this works as preshower
// when PMD is placed in -ve z.
// SUPB --> Pb converter
// SUFE --> Fe backing
// SUPR --> supermodule without rotation (this acts as CPV).
// 
  
  const Int_t npad2 = 72;
  Float_t hexd1[10] = {0.,360.,6,2,-0.4,0.,0.53,0.4,0.,0.53};
//total wall thickness=0.2*2
  Float_t hexd2[10] = {0.,360.,6,2,-0.4,0.,0.51,0.4,0.,0.51};
  Int_t i, j;
  Float_t xb, yb, zb;//, sw[3];
  Int_t number;
  Int_t ihrotm,irotdm;
  const Float_t root3_cons = sqrt(3) /2.; 
  Int_t *idtmed = fIdtmed->GetArray()-599;
 
  AliMatrix(ihrotm, 90., 30.,   90.,  120., 0., 0.);
  AliMatrix(irotdm, 90., 180.,  90.,  270., 180., 0.);
  zdist1  = fIn[2];
  zdist = TMath::Abs(zdist1);
//
  Int_t xrow=1;
  Float_t dpara[6] = {12.5,12.5,0.4,30.,0.,0.};
  dpara[0]=(npad2+0.25)*hexd1[6];
  dpara[1] = dpara[0] *root3_cons;
//
//Subhasis, dimensional parameters of rhombus (dpara) as given to gsvolu
// rhombus to accomodate 72 x 72 hexagons, and with total 1.2cm extension  
//(1mm tolerance on both side and 5mm thick G10 wall)
// 
  
// **** PAD SIZE 10 MM SQUARE EQUIVALENT
//
// Inner hex filled with gas
  gMC->Gsvolu("DP11", "PGON", idtmed[604], hexd2,10);
  gMC->Gsatt("DP11", "SEEN", 1);

// Outer hex filled with Plastic
//plastic  gMC->Gsvolu("DS11", "PGON", idtmed[616], hexd1,10);
// Iron
  gMC->Gsvolu("DS11", "PGON", idtmed[601], hexd1,10);
  gMC->Gsatt("DS11", "SEEN", 1);
// --- place  inner hex inside outer hex 
    gMC->Gsposp("DP11", 1, "DS11", 0., 0., 0., 0, "ONLY", hexd2, 10);
// Rhombus shaped supermodules (defined by PARA)
// volume for SUPERMODULE 
  Float_t dpara_sm[6] = {12.5,12.5,0.8,30.,0.,0.};
  dpara_sm[0]=(npad2+0.25)*hexd1[6] + 1.2;
  dpara_sm[1] = dpara_sm[0] *root3_cons;
//  
  gMC->Gsvolu("SUPR","PARA", idtmed[607], dpara_sm, 6);
  gMC->Gsatt("SUPR", "SEEN", 1);
//  SS 
  Float_t dpara_ss[6] = {12.5,12.5,8.,30.,0.,0.};
  dpara_ss[0]= dpara[0];
  dpara_ss[1]= dpara[1];
  dpara_ss[2]= 0.3/2.;
//
  gMC->Gsvolu("SMSS","PARA", idtmed[601], dpara_ss, 6);
  gMC->Gsatt("SMSS", "SEEN", 1);
// Air 
  Float_t dpara_air[6] = {12.5,12.5,8.,30.,0.,0.};
  dpara_air[0]= dpara[0] - 0.5;
  dpara_air[1]= dpara_air[0] * root3_cons;
  dpara_air[2]= 0.1/2.;
//  gMC->Gsvolu("SMAR","PARA", idtmed[604], dpara_air, 6);
  gMC->Gsvolu("SMAR","PARA", idtmed[698], dpara_air, 6);
  gMC->Gsatt("SMAR", "SEEN", 1);
//  
// volume for gas chamber (DW11)
//  
//  gMC->Gsvolu("DW11","PARA", idtmed[604], dpara, 6);
  gMC->Gsvolu("DW11","PARA", idtmed[698], dpara, 6);
  gMC->Gsatt("DW11", "SEEN", 1);
// Place outer hex inside DW11
  yb = -dpara[1] + (1./root3_cons)*hexd1[6];
  zb = 0.;
  for (j = 1; j <= npad2; ++j) {
  xb =-(dpara[0] + dpara[1]*0.577) + 2*hexd1[6];
   if(xrow >= 2){
    xb = xb+(xrow-1)*hexd1[6];
    }
  for (i = 1; i <= npad2; ++i) {
      number = i+(j-1)*npad2;
    gMC->Gsposp("DS11", number, "DW11", xb, yb, zb, ihrotm, "ONLY", hexd1, 10);
    xb += (hexd1[6]*2.);
  }
   xrow = xrow+1;
    yb += (hexd1[6]*sqrt(3.));
  }
 Float_t z_ss,z_air1,z_air2,z_gas; 
// Place other components inside super module 
    z_ss=-dpara_sm[2]+dpara_ss[2]; 
    gMC->Gspos("SMSS", 1, "SUPR", 0., 0., z_ss, 0, "ONLY");
    z_air1=z_ss+dpara_ss[2] +dpara_air[2]; 
    gMC->Gspos("SMAR", 1, "SUPR", 0., 0., z_air1, 0, "ONLY");
    z_gas=z_air1+dpara_air[2]+dpara[2]+0.1; 
    gMC->Gspos("DW11", 1, "SUPR", 0., 0., z_gas, 0, "ONLY");
    z_air2=z_gas+dpara[2]+0.1+dpara_air[2]; 
    gMC->Gspos("SMAR", 2, "SUPR", 0., 0., z_air2, 0, "ONLY");
  
// --- DEFINE MODules, iron, and lead voLUMES 
  
  
// volume for SUPERMODULE 
//   Pb 
  Float_t dpara_pb[6] = {12.5,12.5,8.,30.,0.,0.};
  dpara_pb[0]=dpara_sm[0];
  dpara_pb[1]=dpara_sm[1];
//  dpara_pb[2]=1.1/2.;
  dpara_pb[2]=1.5/2.;
  gMC->Gsvolu("SUPB","PARA", idtmed[600], dpara_pb, 6);
  gMC->Gsatt("SUPB", "SEEN", 1);
//   Fe 
  Float_t dpara_fe[6] = {12.5,12.5,8.,30.,0.,0.};
  dpara_fe[0]=dpara_sm[0];
  dpara_fe[1]=dpara_sm[1];
  dpara_fe[2]=0.5/2.;
  gMC->Gsvolu("SUFE","PARA", idtmed[601], dpara_fe, 6);
  gMC->Gsatt("SUFE", "SEEN", 1);
// volume for DM11 
  Float_t dpara_dm11[6] = {12.5,12.5,0.8,30.,0.,0.};
  dpara_dm11[0]=dpara_sm[0]+.01;
  dpara_dm11[1] = dpara_dm11[0] *root3_cons;
  dpara_dm11[2]= 6.2/2.;

//  
  gMC->Gsvolu("DM11","PARA", idtmed[698], dpara_dm11, 6);
  gMC->Gsatt("DM11", "SEEN", 1);
// position super module inside DM11
 Float_t z_ps,z_pb,z_fe,z_cv; 
  z_ps=-dpara_dm11[2]+dpara_sm[2];
  gMC->Gspos("SUPR", 1, "DM11", 0., 0., z_ps, irotdm, "ONLY");
  z_pb=z_ps+dpara_sm[2]+dpara_pb[2];
  gMC->Gspos("SUPB", 1, "DM11", 0., 0., z_pb, 0, "ONLY");
  z_fe=z_pb+dpara_pb[2]+dpara_fe[2];
  gMC->Gspos("SUFE", 1, "DM11", 0., 0., z_fe, 0, "ONLY");
  z_cv=z_fe+dpara_fe[2]+dpara_sm[2];
  gMC->Gspos("SUPR", 2, "DM11", 0., 0., z_cv, 0, "ONLY");
// 
}
 
//_____________________________________________________________________________
void AliPMDv1::DrawModule()
{
  //
  // Draw a shaded view of the Photon Multiplicity Detector
  //

  gMC->Gsatt("*", "seen", -1);
  gMC->Gsatt("alic", "seen", 0);
  //
  // Set the visibility of the components
  // 
  gMC->Gsatt("DP11","seen",0);
  gMC->Gsatt("DS11","seen",1);
  gMC->Gsatt("DW11","seen",0);
  gMC->Gsatt("DM11","seen",1);
  gMC->Gsatt("DPMD","seen",0);
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
}

//_____________________________________________________________________________
void AliPMDv1::CreateMaterials()
{
  //
  // Create materials for the PMD version 1
  //
  // ORIGIN    : Y. P. VIYOGI 
  //
  
  // --- The Argon- CO2 mixture --- 
  Float_t ag[2] = { 39.95 };
  Float_t zg[2] = { 18. };
  Float_t wg[2] = { .8,.2 };
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
  
  AliMaterial(96, "MYLAR$", 8.73, 4.55, 1.39, 28.7, 62.);
  AliMaterial(97, "CONCR$", 20., 10., 2.5, 10.7, 40.);
  AliMaterial(98, "Vacum$", 1e-9, 1e-9, 1e-9, 1e16, 1e16);
  AliMaterial(99, "Air  $", 14.61, 7.3, .0012, 30420., 67500.);
  AliMixture(19, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel); 
  // 	define gas-mixtures 
  
  char namate[21];
  gMC->Gfmate((*fIdmate)[3], namate, a, z, d, radl, absl, buf, nbuf);
  ag[1] = a;
  zg[1] = z;
  dg = (dar * 4 + dco) / 5;
  AliMixture(5, "ArCO2$", ag, zg, dg, 2, wg);
  
  // Define tracking media 
  AliMedium(1, "Pb conv.$", 1,  0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
  AliMedium(2, " S steel$", 19, 0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
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
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" PMD_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  printf("                 PMD simulation package (v1) initialised\n");
  printf(" parameters of pmd\n");
  printf("%6d %10.2f %10.2f %10.2f %10.2f %10.2f\n",kdet,thmin,thmax,zdist,thlow,thhigh);
  //
  for(i=0;i<80;i++) printf("*");
  printf("\n");
  //
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
  Int_t   vol[5];
  //  char *namep;
  
  if(gMC->GetMedium() == fMedSens && (destep = gMC->Edep())) {
    
    gMC->CurrentVolID(copy);
//    namep=gMC->CurrentVolName();
//    printf("Current vol is %s \n",namep);
    vol[0]=copy;
    gMC->CurrentVolOffID(1,copy);
//    namep=gMC->CurrentVolOffName(1);
//    printf("Current vol 11 is %s \n",namep);
    vol[1]=copy;
    gMC->CurrentVolOffID(2,copy);
//    namep=gMC->CurrentVolOffName(2);
//    printf("Current vol 22 is %s \n",namep);
    vol[2]=copy;
//	if(strncmp(namep,"DW11",4))vol[2]=1;
    gMC->CurrentVolOffID(3,copy);
//    namep=gMC->CurrentVolOffName(3);
//    printf("Current vol 33 is %s \n",namep);
    vol[3]=copy;
    gMC->CurrentVolOffID(4,copy);
//    namep=gMC->CurrentVolOffName(4);
//    printf("Current vol 44 is %s \n",namep);
    vol[4]=copy;
//	printf("volume number %d,%d,%d,%d,%d \n",vol[0],vol[1],vol[2],vol[3],vol[4]);
    gMC->Gdtom(center,hits,1);
    hits[3] = destep*1e9; //Number in eV
    AddHit(gAlice->CurrentTrack(), vol, hits);
  }
}

  
