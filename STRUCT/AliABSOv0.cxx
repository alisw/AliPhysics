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
Revision 1.1  2000/01/12 15:39:30  morsch
Standar version of ABSO

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Muon ABSOrber                                                            //
//  This class contains the description of the muon absorber geometry        //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliABSOClass.gif">
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
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliABSOv0.h"
#include "AliRun.h"
#include "AliConst.h"
 
ClassImp(AliABSOv0)
 
//_____________________________________________________________________________
AliABSOv0::AliABSOv0()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliABSOv0::AliABSOv0(const char *name, const char *title)
       : AliABSO(name,title)
{
  //
  // Standard constructor
  //
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliABSOv0::CreateGeometry()
{
  //
  // Creation of the geometry of the muon absorber
  //
  //Begin_Html
  /*
    <img src="picts/AliABSOv0Tree.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliABSOv0.gif">
  */
  //End_Html


  Int_t *idtmed = fIdtmed->GetArray()-1599;
  
  Float_t par[24], cpar[5], cpar0[5], pcpar[12], tpar[3], tpar0[3]; 
  Float_t dz;
#include "ShieldConst.h"
// Mother volume and outer shielding: Pb
    
  par[0]  = 0.;
  par[1]  = 360.;
  par[2]  = 7.;

  par[3]  = -(abs_l-abs_d)/2.;
  par[4]  = r_abs;
  par[5]  = abs_d * TMath::Tan(theta1);

  par[6]  = par[3]+(z_nose-abs_d);
  par[7]  = r_abs;
  par[8]  = z_nose * TMath::Tan(theta1);

  par[9]  = par[3]+(z_cone-abs_d);
  par[10] = r_abs;
  par[11] = par[8] + (par[9] - par[6]) * TMath::Tan(theta2);

  par[12]  = par[3]+(abs_c-abs_d);
  par[13] = r_abs;
  par[14] = par[11] + (par[12] - par[9]) * TMath::Tan(acc_max);

  par[15] = par[3]+(abs_l-d_rear-abs_d);
  par[16] = r_abs   + (par[15] - par[12]) * TMath::Tan(theta_open1) ;
  par[17] = par[14] + (par[15] - par[12]) * TMath::Tan(acc_max);

  par[18] = par[3]+(abs_l-d_rear-abs_d);
  par[19] = (abs_l-d_rear) * TMath::Tan(acc_min);
  par[20] = par[14] + (par[15] - par[12]) * TMath::Tan(acc_max);

  par[21] = -par[3];
  par[22] =  abs_l* TMath::Tan(acc_min);
  par[23] = par[20] + (par[21] - par[18]) * TMath::Tan(acc_max);
  gMC->Gsvolu("ABSS", "PCON", idtmed[1612], par, 24);
  for (Int_t i=4; i<24; i+=3) par[i]  = 0;
  gMC->Gsvolu("ABSM", "PCON", idtmed[1655], par, 24);
  gMC->Gspos("ABSS", 1, "ABSM", 0., 0., 0., 0, "ONLY");

//
// Steel envelope
//
  par[4] = par[5] -d_steel;
  par[7] = par[8] -d_steel;
  par[10]= par[11]-d_steel;  
  par[13]= par[14]-d_steel;  
  par[16]= par[17]-d_steel;  
  par[19]= par[20]-d_steel;  
  par[22]= par[23]-d_steel;  
  gMC->Gsvolu("ABST", "PCON", idtmed[1618], par, 24);
  gMC->Gspos("ABST", 1, "ABSS", 0., 0., 0., 0, "ONLY");
//
// Polyethylene shield
// 
  cpar[0] = (abs_l - z_cone) / 2.;
  cpar[1] = z_cone * TMath::Tan(acc_max);
  cpar[2] = cpar[1] + d_poly;
  cpar[3] = abs_l * TMath::Tan(acc_max);
  cpar[4] = cpar[3] + d_poly;
  gMC->Gsvolu("APOL", "CONE", idtmed[1657], cpar, 5);
  dz = (abs_l-abs_d)/2.-cpar[0];
  gMC->Gspos("APOL", 1, "ABSS", 0., 0., dz, 0, "ONLY");

//
// Tungsten nose to protect TPC
// 
  cpar[0] = (z_nose - abs_d) / 2.;
  cpar[1] = abs_d * TMath::Tan(acc_max);
  cpar[2] = abs_d * TMath::Tan(theta1)-d_steel;
  cpar[3] = z_nose * TMath::Tan(acc_max);
  cpar[4] = z_nose * TMath::Tan(theta1)-d_steel;
  gMC->Gsvolu("ANOS", "CONE", idtmed[1611], cpar, 5);
  //
  dz = -(abs_l-abs_d)/2.+cpar[0];
  gMC->Gspos("ANOS", 1, "ABSS", 0., 0., dz, 0, "ONLY");
//
// Tungsten inner shield
//
  cpar[0] = (abs_l-d_rear - abs_c)/ 2.;
  cpar[1] = r_abs;
  cpar[2] = abs_c * TMath::Tan(acc_min);
  cpar[3] = r_abs + 2. * cpar[0] * TMath::Tan(theta_open1);
  cpar[4] = (abs_l-d_rear)  * TMath::Tan(acc_min);
  gMC->Gsvolu("AWIN", "CONE", idtmed[1651], cpar, 5);
  //
  dz = (abs_l-abs_d)/2.-cpar[0]-d_rear;
  gMC->Gspos("AWIN", 1, "ABSS", 0., 0., dz, 0, "ONLY");

  //     Inner tracking region
  //
  //     mother volume: Pb
  //
  pcpar[0]  = 0.;
  pcpar[1]  = 360.;
  pcpar[2]  = 3.;
  pcpar[3]  = -(abs_l-abs_d)/2.;
  pcpar[4]  = r_abs;
  pcpar[5]  = abs_d * TMath::Tan(acc_max);
  pcpar[6]  = pcpar[3]+(z_2deg-abs_d);
  pcpar[7]  = r_abs;
  pcpar[8]  = z_2deg * TMath::Tan(acc_max);
  pcpar[9]  = -par[3];
  pcpar[10] = abs_l * TMath::Tan(acc_min);
  pcpar[11] = abs_l * TMath::Tan(acc_max);
  gMC->Gsvolu("AITR", "PCON", idtmed[1612], pcpar, 12);
  //
  // special Pb medium for last 5 cm of Pb
  zr=abs_l-5;
  cpar[0] = 2.5;
  cpar[1] = zr * TMath::Tan(theta_r);
  cpar[2] = zr * TMath::Tan(acc_max);
  cpar[3] = cpar[1] + TMath::Tan(acc_min) * 5;
  cpar[4] = cpar[2] + TMath::Tan(acc_max) * 5;
  gMC->Gsvolu("ARPB", "CONE", idtmed[1632], cpar, 5);
  dz=(abs_l-abs_d)/2.-cpar[0];
  gMC->Gspos("ARPB", 1, "AITR", 0., 0., dz, 0, "ONLY");

  //
  //     concrete cone: concrete 
  //
  pcpar[9]  = par[3]+(abs_l-d_rear-abs_d);
  pcpar[10] = (abs_l-d_rear) * TMath::Tan(acc_min);
  pcpar[11] = (abs_l-d_rear) * TMath::Tan(acc_max);
  gMC->Gsvolu("ACON", "PCON", idtmed[1616], pcpar, 12);
  gMC->Gspos("ACON", 1, "AITR", 0., 0., 0., 0, "ONLY");
  //
  //     carbon cone: carbon
  //
  pcpar[9]  = pcpar[3]+(abs_cc-abs_d);
  pcpar[10]  = abs_cc * TMath::Tan(acc_min);
  pcpar[11]  = abs_cc * TMath::Tan(acc_max);
  gMC->Gsvolu("ACAR", "PCON", idtmed[1605], pcpar, 12);
  gMC->Gspos("ACAR", 1, "ACON", 0., 0., 0., 0, "ONLY");
  //
  //     inner W shield
  zr=abs_l-d_rear;
  cpar[0] = d_rear/2.;
  cpar[1] = zr * TMath::Tan(acc_min);
  cpar[2] = zr * TMath::Tan(theta_r);
  cpar[3] = cpar[1] + TMath::Tan(acc_min) * 35;
  cpar[4] = cpar[2] + TMath::Tan(theta_r) * 35;
  gMC->Gsvolu("ARW0", "CONE", idtmed[1611], cpar, 5);
  dz=(abs_l-abs_d)/2.-cpar[0];
  gMC->Gspos("ARW0", 1, "AITR", 0., 0., dz, 0, "ONLY");
  //
  // special W medium for last 5 cm of W
  zr=abs_l-5;
  cpar[0] = 2.5;
  cpar[1] = zr * TMath::Tan(acc_min);
  cpar[2] = zr * TMath::Tan(theta_r);
  cpar[3] = cpar[1] + TMath::Tan(acc_min) * 5.;
  cpar[4] = cpar[2] + TMath::Tan(theta_r) * 5.;
  gMC->Gsvolu("ARW1", "CONE", idtmed[1631], cpar, 5);
  dz=d_rear/2.-cpar[0];
  gMC->Gspos("ARW1", 1, "ARW0", 0., 0., dz, 0, "ONLY");
  //
  // PolyEthylene Layers
  Float_t dr_min=TMath::Tan(theta_r) * 5;
  Float_t dr_max=TMath::Tan(acc_max) * 5;
  gMC->Gsvolu("ARPE", "CONE", idtmed[1617], cpar, 0);
  cpar[0]=2.5;
  for (Int_t i=0; i<3; i++) {
      zr=abs_l-d_rear+5+i*10.;
      cpar[1] = zr * TMath::Tan(theta_r);
      cpar[2] = zr * TMath::Tan(acc_max);
      cpar[3] = cpar[1] + dr_min;
      cpar[4] = cpar[2] + dr_max;
      dz=(abs_l-abs_d)/2.-cpar[0]-5.-(2-i)*10;
      gMC->Gsposp("ARPE", i+1, "AITR", 0., 0., dz, 0, "ONLY",cpar,5);
  }
  gMC->Gspos("AITR", 1, "ABSS", 0., 0., 0., 0, "ONLY");	
  dz = (abs_l-abs_d)/2.+abs_d;
  gMC->Gspos("ABSM", 1, "ALIC", 0., 0., dz, 0, "ONLY");	
//
//
// vacuum system
//
// pipe and heating jackets
//
//
// cylindrical piece
  tpar0[2]=(abs_c-abs_d)/2;
  tpar0[0]=r_vacu;
  tpar0[1]=r_abs;
  gMC->Gsvolu("AV11", "TUBE", idtmed[1658], tpar0, 3);
//
// insulation
  tpar[2]=tpar0[2];
  tpar[0]=tpar0[0]+d_tube;
  tpar[1]=tpar0[0]+d_tube+d_insu;
  gMC->Gsvolu("AI11", "TUBE", idtmed[1653], tpar, 3);
  gMC->Gspos("AI11", 1, "AV11", 0., 0., 0., 0, "ONLY"); 
//
// clearance
  tpar[0]=tpar0[1]-d_prot-d_free;
  tpar[1]=tpar0[1]-d_prot;
  gMC->Gsvolu("AP11", "TUBE", idtmed[1655], tpar, 3);
  gMC->Gspos("AP11", 1, "AV11", 0., 0., 0., 0, "ONLY"); 
  
  dz=-(abs_l-abs_d)/2.+tpar0[2];
  gMC->Gspos("AV11", 1, "ABSM", 0., 0., dz, 0, "ONLY"); 
  
//
// conical piece
  cpar0[0]=(abs_l-d_rear-abs_c)/2;
  cpar0[1]=r_vacu;
  cpar0[2]=r_abs;
  cpar0[3]=cpar0[1]+2.*cpar0[0]*TMath::Tan(theta_open1);
  cpar0[4]=cpar0[2]+2.*cpar0[0]*TMath::Tan(theta_open1);
  gMC->Gsvolu("AV21", "CONE", idtmed[1658], cpar0, 5);
//
// insulation
  cpar[0]=cpar0[0];
  cpar[1]=cpar0[1]+d_tube;
  cpar[2]=cpar0[1]+d_tube+d_insu;
  cpar[3]=cpar0[3]+d_tube;
  cpar[4]=cpar0[3]+d_tube+d_insu;
  gMC->Gsvolu("AI21", "CONE", idtmed[1653], cpar, 5);
  gMC->Gspos("AI21", 1, "AV21", 0., 0., 0., 0, "ONLY"); 
//
// clearance
  cpar[1]=cpar0[2]-d_prot-d_free;
  cpar[2]=cpar0[2]-d_prot;
  cpar[3]=cpar0[4]-d_prot-d_free;
  cpar[4]=cpar0[4]-d_prot;
  gMC->Gsvolu("AP21", "CONE", idtmed[1655], cpar, 5);
  gMC->Gspos("AP21", 1, "AV21", 0., 0., 0., 0, "ONLY"); 
  
  dz=(abs_l-abs_d)/2.-cpar0[0]-d_rear;
  gMC->Gspos("AV21", 1, "ABSM", 0., 0., dz, 0, "ONLY"); 

}

//_____________________________________________________________________________

void AliABSOv0::Init()
{
  //
  // Initialisation of the muon absorber after it has been built
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" ABSOv0_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}
 









