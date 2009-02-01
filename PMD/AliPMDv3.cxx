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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Photon Multiplicity Detector Version 1                                   //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliPMDv3Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
////

#include "AliPMDv3.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliMC.h" 
#include "AliConst.h" 
#include "iostream.h"
 
static Int_t   kdet, ncell_sm, ncell_hole;
static Float_t zdist, zdist1;
static Float_t sm_length, sm_thick, cell_radius, cell_wall, cell_depth;
static Float_t boundary, th_base, th_air, th_pcb;
static Float_t th_lead, th_steel;

ClassImp(AliPMDv3)
 
  //_____________________________________________________________________________
  AliPMDv3::AliPMDv3()
{
  //
  // Default constructor 
  //
  fMedSens=0;
}
 
//_____________________________________________________________________________
AliPMDv3::AliPMDv3(const char *name, const char *title)
  : AliPMD(name,title)
{
  //
  // Standard constructor
  //
  fMedSens=0;
}

//_____________________________________________________________________________
void AliPMDv3::CreateGeometry()
{
  //
  // Create geometry for Photon Multiplicity Detector Version 3 :
  // April 2, 2001
  //
  //Begin_Html
  /*
    <img src="picts/AliPMDv3.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliPMDv3Tree.gif">
  */
  //End_Html
  GetParameters();
  CreateSupermodule();
  CreatePMD();
}

//_____________________________________________________________________________
void AliPMDv3::CreateSupermodule()
{
  //
  // Creates the geometry of the cells, places them in  supermodule which
  // is a rhombus object.

  // *** DEFINITION OF THE GEOMETRY OF THE PMD  *** 
  // *** HEXAGONAL CELLS WITH 10 MM SQUARE EQUIVALENT
  // -- Author :     S. Chattopadhyay, 02/04/1999. 

  // Basic unit is ECAR, a hexagonal cell made of Ar+CO2, which is placed inside another 
  // hexagonal cell made of Cu (ECCU) with larger radius, compared to ECAR. The difference
  // in radius gives the dimension of half width of each cell wall.
  // These cells are placed as 72 x 72 array in a 
  // rhombus shaped supermodule (EHC1). The rhombus shaped modules are designed
  // to have closed packed structure.
  //
  // Each supermodule (ESM1 or ESM2), made of G10 is filled with following components
  //  EAIR --> Air gap between gas hexagonal cells and G10 backing.
  //  EHC1 --> Rhombus shaped parallelopiped containing the hexagonal cells
  //  EAIR --> Air gap between gas hexagonal cells and G10 backing.
  //
  // ESM1 is placed in EMM1 along with EMPB (Pb converter) and EMFE (iron support) 
  //  EMM1 made of
  // ESM1 --> Normal supermodule
  // EMPB --> Pb converter
  // EMFE --> Fe backing
  //
  // ESM2 is placed in EMM2 along with EMPB (Pb converter) and EMFE (iron support) 
  //  EMM2 made of
  // ESM2 --> Special supermodule containing the cut for the hole
  // EMPB --> Pb converter
  // EMFE --> Fe backing
  
  //
  //  		                     EPMD
  //  				       |             
  //   				       |
  //        -------------------------------------------------------------------
  //        |                   |                          |                   |
  //       EHOL                EMM1                       EMM2               EALM
  //                            |                           |
  //                ----------------------       ------------------------
  //                |       |      |     |       |      |        |      |
  //               ESM1    EMPB  EMFE  ESM1     ESM2   EMPB     EMFE   ESM2
  //                |                            |
  //          ------------                 -------------
  //          |     |     |                |     |      |
  //        EAIR  EHC1   EAIR             EAIR  EHC2  EAIR
  //                |                            |
  //              ECCU                          ECCU
  //                |                            |
  //              ECAR                          ECAR 
  

  Int_t i, j;
  Float_t xb, yb, zb;
  Int_t number;
  Int_t ihrotm,irotdm;
  const Float_t root3_2 = TMath::Sqrt(3.) /2.; 
  Int_t *idtmed = fIdtmed->GetArray()-599;
 
  AliMatrix(ihrotm, 90., 30.,   90.,  120., 0., 0.);
  AliMatrix(irotdm, 90., 180.,  90.,  270., 180., 0.);
 
  zdist = TMath::Abs(zdist1);


  //Subhasis, dimensional parameters of rhombus (dpara) as given to gsvolu
  // rhombus to accomodate 72 x 72 hexagons, and with total 1.2cm extension  
  //(1mm tolerance on both side and 5mm thick G10 wall)
  // 
  
  // **** CELL SIZE 20mm^2 EQUIVALENT

  // Inner hexagon filled with gas (Ar+CO2)

  Float_t hexd2[10] = {0.,360.,6,2,-0.25,0.,0.23,0.25,0.,0.23};

  hexd2[4]= - cell_depth/2.;
  hexd2[7]=   cell_depth/2.;
  hexd2[6]=   cell_radius - cell_wall;
  hexd2[9]=   cell_radius - cell_wall;
  
  gMC->Gsvolu("ECAR", "PGON", idtmed[604], hexd2,10);
  gMC->Gsatt("ECAR", "SEEN", 0);

  // Outer hexagon made of Copper

  Float_t hexd1[10] = {0.,360.,6,2,-0.25,0.,0.25,0.25,0.,0.25};
  //total wall thickness=0.2*2

  hexd1[4]= - cell_depth/2.;
  hexd1[7]=   cell_depth/2.;
  hexd1[6]=   cell_radius;
  hexd1[9]=   cell_radius;

  gMC->Gsvolu("ECCU", "PGON", idtmed[614], hexd1,10);
  gMC->Gsatt("ECCU", "SEEN", 1);


// Rhombus shaped supermodules (defined by PARA) 

// volume for SUPERMODULE 
   
  Float_t dpara_sm[6] = {12.5,12.5,0.8,30.,0.,0.};
  dpara_sm[0]=(ncell_sm+0.25)*hexd1[6] ;
  dpara_sm[1] = dpara_sm[0] *root3_2;
  dpara_sm[2] = sm_thick/2.;

//  G10 inner part of supermodule, these will be 9 in all, one being special

  Float_t dpara_g10[6] = {12.5,12.5,8.,30.,0.,0.};
  dpara_g10[0]= dpara_sm[0];
  dpara_g10[1]= dpara_sm[1];
  dpara_g10[2]= dpara_sm[2];

//
  gMC->Gsvolu("ESM1","PARA", idtmed[607], dpara_g10, 6);
  gMC->Gsatt("ESM1", "SEEN", 0);
  //
  gMC->Gsvolu("ESM2","PARA", idtmed[607], dpara_g10, 6);
  gMC->Gsatt("ESM2", "SEEN", 0);

  // Air residing between the PCB and the base

  Float_t dpara_air[6] = {12.5,12.5,8.,30.,0.,0.};
  dpara_air[0]= dpara_sm[0];
  dpara_air[1]= dpara_sm[1];
  dpara_air[2]= th_air/2.;

  gMC->Gsvolu("EAIR","PARA", idtmed[698], dpara_air, 6);
  gMC->Gsatt("EAIR", "SEEN", 0);

  // volume for honeycomb chamber (EHC1 and EHC2)

  Float_t dpara[6] = {12.5,12.5,0.4,30.,0.,0.};
  dpara[0] = dpara_sm[0];
  dpara[1] = dpara_sm[1];
  dpara[2] = cell_depth/2.;

  gMC->Gsvolu("EHC1","PARA", idtmed[698], dpara, 6);
  gMC->Gsatt("EHC1", "SEEN", 1);
  
  gMC->Gsvolu("EHC2","PARA", idtmed[698], dpara, 6);
  gMC->Gsatt("EHC2", "SEEN", 1);

  // --- place  inner hex inside outer hex 

  gMC->Gsposp("ECAR", 1, "ECCU", 0., 0., 0., 0, "ONLY", hexd2, 10);

  // Place outer hex ECCU cells  inside EHC1 (72 X 72)

  Int_t xrow=1;

  yb = -dpara[1] + (1./root3_2)*hexd1[6];
  zb = 0.;

  for (j = 1; j <= ncell_sm; ++j) {
    xb =-(dpara[0] + dpara[1]*0.577) + 2*hexd1[6]; //0.577=tan(30deg)
    if(xrow >= 2){
      xb = xb+(xrow-1)*hexd1[6];
    }
    for (i = 1; i <= ncell_sm; ++i) {
      number = i+(j-1)*ncell_sm;
      gMC->Gsposp("ECCU", number, "EHC1", xb,yb,zb, ihrotm, "ONLY", hexd1,10);
      xb += (hexd1[6]*2.);
    }
    xrow = xrow+1;
    yb += (hexd1[6]*TMath::Sqrt(3.));
  }


  // Place outer hex ECCU inside EHC2
  // skip cells which go into the hole in top left corner.

  xrow=1;
  yb = -dpara[1] + (1./root3_2)*hexd1[6];
  zb = 0.;
  for (j = 1; j <= ncell_sm; ++j) {
    xb =-(dpara[0] + dpara[1]*0.577) + 2*hexd1[6];
    if(xrow >= 2){
      xb = xb+(xrow-1)*hexd1[6];
    }
    for (i = 1; i <= ncell_sm; ++i) {
      number = i+(j-1)*ncell_sm;
      if(i > ncell_hole || j <= (ncell_sm - ncell_hole))
	{
	  gMC->Gsposp("ECCU", number, "EHC2", xb,yb,zb, ihrotm, "ONLY", hexd1,10);
	} 
      xb += (hexd1[6]*2.);
    }
    xrow = xrow+1;
    yb += (hexd1[6]*TMath::Sqrt(3.));
  }

  // Place EHC1 and EAIR into  ESM1; EHC2 and EAIR into ESM2 

  Float_t z_air1,z_air2,z_gas; 

  z_air1= -dpara_g10[2] + th_base + dpara_air[2]; 
  gMC->Gspos("EAIR", 1, "ESM1", 0., 0., z_air1, 0, "ONLY");
  z_gas=z_air1+dpara_air[2]+ th_pcb + dpara[2]; 
  gMC->Gspos("EHC1", 1, "ESM1", 0., 0., z_gas, 0, "ONLY");
  z_air2=z_gas+dpara[2]+ th_pcb + dpara_air[2]; 
  gMC->Gspos("EAIR", 2, "ESM1", 0., 0., z_air2, 0, "ONLY");

  z_air1= -dpara_g10[2] + th_base + dpara_air[2]; 
  gMC->Gspos("EAIR", 1, "ESM2", 0., 0., z_air1, 0, "ONLY");
  z_gas=z_air1+dpara_air[2]+ th_pcb + dpara[2]; 
  gMC->Gspos("EHC2", 1, "ESM2", 0., 0., z_gas, 0, "ONLY");
  z_air2=z_gas+dpara[2]+ th_pcb + dpara_air[2]; 
  gMC->Gspos("EAIR", 2, "ESM2", 0., 0., z_air2, 0, "ONLY");

}
 
//_____________________________________________________________________________

void AliPMDv3::CreatePMD()
{
  //
  // Create final detector from supermodules
  //
  // -- Author :     Y.P. VIYOGI, 07/05/1996. 
  // -- Modified:    P.V.K.S.Baba(JU), 15-12-97. 
  // -- Modified:    For New Geometry YPV, March 2001.


  // Gaspmd, the dimension of TUBE mother volume of PMD,

  Float_t gaspmd[3] = { 0.,150.,10.};

  const Float_t root3_2 = TMath::Sqrt(3.)/2.;
  const Float_t pi = 3.14159;
  Int_t i,j;

  Float_t  xp, yp, zp;

  Int_t num_mod;
  Int_t jhrot12,jhrot13, irotdm;

  Int_t *idtmed = fIdtmed->GetArray()-599;
  
  //  VOLUMES Names : begining with D for all PMD volumes, 
  // The names of SIZE variables begin with S and have more meaningful
  // characters as shown below. 
  
  // 		VOLUME 	SIZE	MEDIUM	: 	REMARKS 
  // 		------	-----	------	: --------------------------- 
  
  // 		EPMD	GASPMD	 AIR	: INSIDE PMD  and its SIZE 
  
  // *** Define the  EPMD   Volume and fill with air *** 

  gMC->Gsvolu("EPMD", "TUBE", idtmed[698], gaspmd, 3);
  gMC->Gsatt("EPMD", "SEEN", 0);

  AliMatrix(irotdm, 90., 0.,  90.,  90., 180., 0.);
   
  AliMatrix(jhrot12, 90., 120., 90., 210., 0., 0.);
  AliMatrix(jhrot13, 90., 240., 90., 330., 0., 0.);

  // dpara_emm1 array contains parameters of the imaginary volume EMM1, 
  // this is just a little more than the side of a supermodule. 

  Float_t dm_thick = 2. * sm_thick + th_lead + th_steel;

  Float_t dpara_emm1[6] = {12.5,12.5,0.8,30.,0.,0.};
  dpara_emm1[0] = sm_length/2.;
  dpara_emm1[1] = dpara_emm1[0] *root3_2;
  dpara_emm1[2] = dm_thick/2.;

  // EMM1 : normal volume as in old cases
  gMC->Gsvolu("EMM1","PARA", idtmed[698], dpara_emm1, 6);
  gMC->Gsatt("EMM1", "SEEN", 1);

  // EMM2 : special volume containing special supermodule
  gMC->Gsvolu("EMM2","PARA", idtmed[698], dpara_emm1, 6);
  gMC->Gsatt("EMM2", "SEEN", 1);

  //
  // --- DEFINE MODules, iron, and lead voLUMES 
  
  //place ESM1 into EMM1 and ESM2 into EMM2 along with EMPB and EMFE

  Float_t dx = sm_length;
  Float_t dy = dx * root3_2;

  Float_t xsup[9] = {-dx/2., dx/2., 3.*dx/2., 
		     -dx,    0.,       dx,
		     -3.*dx/2., -dx/2., dx/2.};

  Float_t ysup[9] = {dy,  dy,  dy, 
                     0.,  0.,  0., 
                    -dy, -dy, -dy};

  //
  
  // volume for SUPERMODULE 

  //   Pb Convertor 
  Float_t dpara_pb[6] = {12.5,12.5,8.,30.,0.,0.};
  dpara_pb[0] = sm_length/2.;
  dpara_pb[1] = dpara_pb[0] * root3_2;
  dpara_pb[2] = th_lead/2.;

  gMC->Gsvolu("EMPB","PARA", idtmed[600], dpara_pb, 6);
  gMC->Gsatt ("EMPB", "SEEN", 0);

  //   Fe Support
  Float_t dpara_fe[6] = {12.5,12.5,8.,30.,0.,0.};
  dpara_fe[0] = dpara_pb[0];
  dpara_fe[1] = dpara_pb[1];
  dpara_fe[2] = th_steel/2.;

  gMC->Gsvolu("EMFE","PARA", idtmed[618], dpara_fe, 6);
  gMC->Gsatt ("EMFE", "SEEN", 0);

  // position supermodule ESM1 inside EMM1

  Float_t z_ps,z_pb,z_fe,z_cv; 
  
  z_ps = - dpara_emm1[2] + sm_thick/2.;
  gMC->Gspos("ESM1", 2, "EMM1", 0., 0., z_ps, irotdm, "ONLY");
  z_pb=z_ps+sm_thick/2.+dpara_pb[2];
  gMC->Gspos("EMPB", 1, "EMM1", 0., 0., z_pb, 0, "ONLY");
  z_fe=z_pb+dpara_pb[2]+dpara_fe[2];
  gMC->Gspos("EMFE", 1, "EMM1", 0., 0., z_fe, 0, "ONLY");
  z_cv=z_fe+dpara_fe[2]+sm_thick/2.;
  gMC->Gspos("ESM1", 1, "EMM1", 0., 0., z_cv, 0, "ONLY");

  // position supermodule  ESM2 inside EMM2

  z_ps = - dpara_emm1[2] + sm_thick/2.;
  gMC->Gspos("ESM2", 2, "EMM2", 0., 0., z_ps, irotdm, "ONLY");
  z_pb = z_ps + sm_thick/2.+dpara_pb[2];
  gMC->Gspos("EMPB", 1, "EMM2", 0., 0., z_pb, 0, "ONLY");
  z_fe = z_pb + dpara_pb[2]+dpara_fe[2];
  gMC->Gspos("EMFE", 1, "EMM2", 0., 0., z_fe, 0, "ONLY");
  z_cv = z_fe + dpara_fe[2]+sm_thick/2.;
  gMC->Gspos("ESM2", 1, "EMM2", 0., 0., z_cv, 0, "ONLY");
  // 

  // EHOL is a tube structure made of air

  Float_t d_hole[3];
  d_hole[0] = 0.;
  d_hole[1] = ncell_hole * cell_radius *2. * root3_2 + boundary;
  d_hole[2] = dm_thick/2.;

  gMC->Gsvolu("EHOL", "TUBE", idtmed[698], d_hole, 3);
  gMC->Gsatt("EHOL", "SEEN", 1);

  //Al-rod as boundary of the supermodules

  Float_t Al_rod[3] ;
  Al_rod[0] = sm_length * 3/2.;
  Al_rod[1] = boundary;
  Al_rod[2] = dm_thick/2.;

  gMC->Gsvolu("EALM","BOX ", idtmed[698], Al_rod, 3);
  gMC->Gsatt ("EALM", "SEEN", 1);
  Float_t xalm[3];
  xalm[0]=Al_rod[0];
  xalm[1]=-xalm[0]/2.;
  xalm[2]=xalm[1];

  Float_t yalm[3];
  yalm[0]=0.;
  yalm[1]=xalm[0]*root3_2;
  yalm[2]=-yalm[1];

  // delx = full side of the supermodule
  Float_t delx=sm_length * 3.;
  Float_t x1= delx*root3_2 /2.;
  Float_t x4=delx/4.; 

  // xpos and ypos are the x & y coordinates of the centres of EMM1 volumes

  Float_t xoff = boundary * TMath::Tan(pi/6.);
  Float_t xmod[3]={x4 + xoff , x4 + xoff, -2.*x4-boundary/root3_2};
  Float_t ymod[3] = {-x1 - boundary, x1 + boundary, 0.};
  Float_t xpos[9], ypos[9];
  Float_t theta[3] = {0., 2.*pi/3., 4.*pi/3.};
  Int_t irotate[3] = {0, jhrot12, jhrot13};

  for (j=0; j<3; ++j)
    {
      gMC->Gsposp("EALM", j+1, "EPMD", xalm[j],yalm[j], 0., irotate[j], "ONLY", Al_rod, 3);
      for (i=0; i<9; ++i)
	{
	  xpos[i]=xmod[j] + xsup[i]*TMath::Cos(theta[j]) - ysup[i]*TMath::Sin(theta[j]);
	  ypos[i]=ymod[j] + xsup[i]*TMath::Sin(theta[j]) + ysup[i]*TMath::Cos(theta[j]);

	  if(fDebug) printf("%s: %f %f \n", ClassName(), xpos[i], ypos[i]);

	  num_mod = i + 1 + j*9;

	  printf("\n%s: Num_mod %d\n",ClassName(),num_mod);

	  if(i==0){ 
	    gMC->Gsposp("EMM2", num_mod, "EPMD", xpos[i],ypos[i], 0., irotate[j], "ONLY", dpara_emm1, 6);
	  }
	  else {
	    gMC->Gsposp("EMM1", num_mod, "EPMD", xpos[i],ypos[i], 0., irotate[j], "ONLY", dpara_emm1, 6);
	  }
	}
    }

	
  // place EHOL in the centre of EPMD
  gMC->Gspos("EHOL", 1, "EPMD", 0.,0.,0., 0, "ONLY");

  // --- Place the EPMD in ALICE 
  xp = 0.;
  yp = 0.;
  zp = zdist1;
  
  gMC->Gspos("EPMD", 1, "ALIC", xp,yp,zp, 0, "ONLY");
    
}

 
//_____________________________________________________________________________
void AliPMDv3::DrawModule()
{
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
  gMC->Gsatt("EHC1","seen",1);
  gMC->Gsatt("EHC1","seen",1);
  gMC->Gsatt("EHC2","seen",1);
  gMC->Gsatt("EMM1","seen",1);
  gMC->Gsatt("EHOL","seen",1);
  gMC->Gsatt("EPMD","seen",0);
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
void AliPMDv3::CreateMaterials()
{
  //
  // Create materials for the PMD
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
  //Float_t x0xe=2.4;
  //Float_t dxe=0.005858;
  Float_t buf[1];
  Int_t nbuf;
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  
  Int_t *idtmed = fIdtmed->GetArray()-599;
  Int_t isxfld = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
  Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();
  
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
  
  char namate[21]="";
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
}

//_____________________________________________________________________________
void AliPMDv3::Init()
{
  //
  // Initialises PMD detector after it has been built
  //
  Int_t i;
  kdet=1;
  //
  if(fDebug) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" PMD_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n");
    printf("%s:                  PMD simulation package (v3) initialised\n",
	   ClassName());
    printf("%s:  parameters of pmd\n",ClassName());
    printf("%s: %10.2f %10.2f %10.2f %10.2f\n",ClassName(),
	   cell_radius,cell_wall,cell_depth,zdist1 );
    printf("%s: ",ClassName());
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
  
  Int_t *idtmed = fIdtmed->GetArray()-599;
  fMedSens=idtmed[605-1];
}

//_____________________________________________________________________________
void AliPMDv3::StepManager()
{
  //
  // Called at each step in the PMD
  //
  Int_t   copy;
  Float_t hits[4], destep;
  Float_t center[3] = {0,0,0};
  Int_t   vol[5];
  //char *namep;
  
  if(gMC->GetMedium() == fMedSens && (destep = gMC->Edep())) {
  
    gMC->CurrentVolID(copy);

    //namep=gMC->CurrentVolName();
    //printf("Current vol is %s \n",namep);

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
    //printf("volume number %d,%d,%d,%d,%d,%f \n",vol[0],vol[1],vol[2],vol[3],vol[4],destep*1000000);

    gMC->Gdtom(center,hits,1);
    hits[3] = destep*1e9; //Number in eV
    AddHit(gAlice->CurrentTrack(), vol, hits);
  }
}

  
//------------------------------------------------------------------------
// Get parameters

void AliPMDv3::GetParameters()
{
  Int_t ncell_um, num_um;
  ncell_um=24;
  num_um=3;
  ncell_hole=12;
  cell_radius=0.25;
  cell_wall=0.02;
  cell_depth=0.25 * 2.;
  //
  boundary=0.7;
  ncell_sm=ncell_um * num_um;  //no. of cells in a row in one supermodule
  sm_length= ((ncell_sm + 0.25 ) * cell_radius) * 2.;
  //
  th_base=0.3;
  th_air=0.1;
  th_pcb=0.16;
  //
  sm_thick = th_base + th_air + th_pcb + cell_depth + th_pcb + th_air + th_pcb;
  //
  th_lead=1.5;
  th_steel=0.5;
  //
  zdist1 = -370.;
}













