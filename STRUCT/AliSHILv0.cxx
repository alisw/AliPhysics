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
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Muon Shield Class                                                        //
//  This class contains a description of the muon shield                     //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliSHILClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliSHILv0.h"
#include "AliRun.h"
#include "AliConst.h"

ClassImp(AliSHILv0)
 
//_____________________________________________________________________________
AliSHILv0::AliSHILv0()
{
  //
  // Default constructor for muon shield
  //
}
 
//_____________________________________________________________________________
AliSHILv0::AliSHILv0(const char *name, const char *title)
  : AliSHIL(name,title)
{
  //
  // Standard constructor for muon shield
  //
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliSHILv0::CreateGeometry()
{
  //
  // Build muon shield geometry
  //
  //
  //Begin_Html
  /*
    <img src="picts/AliSHILv0.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliSHILv0Tree.gif">
  */
  //End_Html

  Float_t cpar[5], cpar0[5], tpar[3], par1[39], par2[21], par3[27], par4[21], par0[42];
  Float_t dz, DZ;
  
  Int_t *idtmed = fIdtmed->GetArray()-1699;

#include "ShieldConst.h"
//
// Mother volume
//
  par0[0]  = 0.;
  par0[1]  = 360.;
  par0[2]  = 13.;

  Float_t dl=(zvac12-abs_l+d_rear)/2.;
  dz=abs_l-d_rear+dl;
//
  par0[3]  = -dl;
  par0[4]  = 0.;
  par0[5]  = (abs_l-d_rear) * TMath::Tan(acc_min);

  par0[6]  = -dl+d_rear;
  par0[7]  = 0.;
  par0[8]  = abs_l * TMath::Tan(acc_min);

  par0[9]  = -dl+d_rear;
  par0[10]  = 0.;
  par0[11]  = R11;

  par0[12]  = -dz+zvac4;
  par0[13]  = 0.;
  par0[14]  = R11;

  par0[15]  = -dz+zvac4;
  par0[16] = 0.;
  par0[17] = R21;

  par0[18] = -dz+zvac6;
  par0[19] = 0.;
  par0[20] = R21;

  par0[21] = -dz+zvac6;
  par0[22] = 0.;
  par0[23] = zvac6 * TMath::Tan(acc_min);

  par0[24] = -dz+zcone_e;
  par0[25] = 0.;
  par0[26] = 30.;

  par0[27] = -dz+zvac10;
  par0[28] = 0.;
  par0[29] = 30.;

  par0[30] = -dz+zvac10;
  par0[31] = 0.;
  par0[32] = R42;

  par0[33] = -dz+zvac11;
  par0[34] = 0.;
  par0[35] = R42;

  par0[36] = -dz+zvac11;
  par0[37] = 0.;
  par0[38] = R43;

  par0[39] = -dz+zvac12;
  par0[40] = 0.;
  par0[41] = R43;

  gMC->Gsvolu("YMOT", "PCON", idtmed[1755], par0, 42);
  dz=abs_l+dl-d_rear;
  gMC->Gspos("YMOT", 1, "ALIC", 0., 0., dz, 0, "ONLY");  
//

  DZ=-dl;

//
// First section: bellows below and behind front absorber 
// 
//

  par1[0]  = 0.;
  par1[1]  = 360.;
  par1[2]  = 12.;
  dl=(zvac4-abs_l+d_rear)/2.;
  
  par1[3]  = -dl;
  par1[4]  = r_abs+(abs_l-d_rear-abs_c) * TMath::Tan(theta_open1);
  par1[5]  = (abs_l-d_rear)* TMath::Tan(acc_min);

  par1[6]  = -dl+zvac1-(abs_l-d_rear);
  par1[7]  = r_abs+ (zvac1-abs_c) * TMath::Tan(theta_open1);
  par1[8]  = zvac1 * TMath::Tan(acc_min);

  par1[9]  = par1[6]+dr11;
  par1[10] = par1[7]+dr11;
  par1[11] = (zvac1+dr11) * TMath::Tan(acc_min);

  par1[12] = -dl+d_rear;
  par1[13] = par1[10];
  par1[14] = abs_l * TMath::Tan(acc_min);

  par1[15] = -dl+d_rear;
  par1[16] = par1[10];
  par1[17] = R11;

  par1[18] = -dl+(zvac1+dr11+dB1-abs_l+d_rear);
  par1[19] = par1[16];
  par1[20] = R11;

  par1[21] = par1[18]+dr12;
  par1[22] = par1[19]+dr12;
  par1[23] = R11;

  par1[24] = par1[21]+dF1;
  par1[25] = par1[22];
  par1[26] = R11;

  par1[27] = par1[24]+dr12;
  par1[28] = par1[25]-dr12;
  par1[29] = R11;

  par1[30] = par1[27]+dB1;
  par1[31] = par1[28];
  par1[32] = R11;

  par1[33] = par1[30]+dr13;
  par1[34] = par1[31]-dr13;
  par1[35] = R11;

  par1[36] = -dl+zvac4-abs_l+d_rear;
  par1[37] = par1[34]+(zvac4-zvac3)*TMath::Tan(theta_open2);
  par1[38] = R11;

  Float_t r2=par1[34];
  
  gMC->Gsvolu("YGO1", "PCON", idtmed[1760], par1, 39);

  for (Int_t i=4; i<38; i+=3) par1[i]  = 0;
  gMC->Gsvolu("YMO1", "PCON", idtmed[1755], par1, 39);

  gMC->Gspos("YGO1", 1, "YMO1", 0., 0., 0., 0, "ONLY");  

  DZ+=dl;
  gMC->Gspos("YMO1", 1, "YMOT", 0., 0., DZ, 0, "ONLY");  
  DZ+=dl;

//
// Steel envelope
  tpar[0]=R11-dRSteel1;
  tpar[1]=R11;
  tpar[2]=dl-d_rear/2;
  gMC->Gsvolu("YSE1", "TUBE", idtmed[1718], tpar, 3);
  dz=dl-tpar[2];
  
  gMC->Gspos("YSE1", 1, "YGO1", 0., 0., dz, 0, "ONLY");  

//
// 1st section: vacuum system
//
//
// Bellow 1
//
  tpar[0]=rB1;
  tpar[1]=rB1+hB1;
  tpar[2]=eB1/2.;
  gMC->Gsvolu("YB11", "TUBE", idtmed[1758], tpar, 3);
  Float_t dl1=tpar[2];
  
  tpar[0]=rB1+hB1-eB1;
  tpar[1]=rB1+hB1;
  tpar[2]=(lB1/2.-2.*eB1)/2.;
  gMC->Gsvolu("YB12", "TUBE", idtmed[1758], tpar, 3);
  Float_t dl2=tpar[2];

  tpar[0]=rB1-eB1;
  tpar[1]=rB1;
  tpar[2]=lB1/8.;
  gMC->Gsvolu("YB13", "TUBE", idtmed[1758], tpar, 3);
  Float_t dl3=tpar[2];


  tpar[0]=0;
  tpar[1]=rB1+hB1;
  tpar[2]=lB1/2.;
  gMC->Gsvolu("YBU1", "TUBE", idtmed[1755], tpar, 3);

  dz=-tpar[2]+dl3;
  gMC->Gspos("YB13", 1, "YBU1", 0., 0., dz, 0, "ONLY"); 
  dz+=dl3;
  dz+=dl1;  
  gMC->Gspos("YB11", 1, "YBU1", 0., 0., dz, 0, "ONLY"); 
  dz+=dl1;  
  dz+=dl2;  
  gMC->Gspos("YB12", 1, "YBU1", 0., 0., dz, 0, "ONLY"); 
  dz+=dl2;  
  dz+=dl1;
  gMC->Gspos("YB11", 2, "YBU1", 0., 0., dz, 0, "ONLY"); 
  dz+=dl1;
  dz+=dl3;
  gMC->Gspos("YB13", 2, "YBU1", 0., 0., dz, 0, "ONLY"); 
  

  tpar[0]=0;
  tpar[1]=rB1+hB1;
  tpar[2]=10.*lB1/2.;
  gMC->Gsvolu("YBM1", "TUBE", idtmed[1755], tpar, 3);
  dz=-tpar[2]+lB1/2.;
  for (Int_t i=0; i<10; i++) {
      gMC->Gspos("YBU1", i+1 , "YBM1", 0., 0.,dz , 0, "ONLY"); 
      dz+=lB1;
  }
  dz=-dl+(zvac1-abs_l+d_rear)+dr11+tpar[2];
  gMC->Gspos("YBM1", 1, "YMO1", 0., 0., dz, 0, "ONLY"); 

  dz=dl-dr13-(zvac4-zvac3)-tpar[2]-(dr11-dr13);
  gMC->Gspos("YBM1", 2, "YMO1", 0., 0., dz, 0, "ONLY"); 

//
// Flange

  tpar[0]=0;
  tpar[1]=rF1;
  tpar[2]=dF1/2.;
  gMC->Gsvolu("YFM1", "TUBE", idtmed[1755], tpar, 3);

  tpar[0]=rF1-d_flange;
  tpar[1]=rF1;
  tpar[2]=dF1/2.;
  gMC->Gsvolu("YF11", "TUBE", idtmed[1758], tpar, 3);
  gMC->Gspos("YF11", 1, "YFM1", 0., 0., 0., 0, "ONLY"); 

  tpar[0]=rB1;
  tpar[1]=rF1-d_flange;
  tpar[2]=d_flange/2.;
  gMC->Gsvolu("YF12", "TUBE", idtmed[1758], tpar, 3);
  dz=-dF1/2.+tpar[2];
  gMC->Gspos("YF12", 1, "YFM1", 0., 0., dz, 0, "ONLY"); 
  dz= dF1/2.-tpar[2];
  gMC->Gspos("YF12", 2, "YFM1", 0., 0., dz, 0, "ONLY"); 

  dz=-dl+(zvac2-abs_l+d_rear);
  gMC->Gspos("YFM1", 2, "YMO1", 0., 0., dz, 0, "ONLY"); 


//
// pipe between flange and bellows
  tpar[0]=rB1-d_tube;
  tpar[1]=rB1;
  tpar[2]=2.*(dB1+dr12-10.*lB1)/4.;
  gMC->Gsvolu("YPF1", "TUBE", idtmed[1758], tpar, 3);
  dz=-dl+(zvac2-abs_l+d_rear)-dF1/2.-tpar[2];
  gMC->Gspos("YPF1", 1, "YMO1", 0., 0., dz, 0, "ONLY"); 
  dz=-dl+(zvac2-abs_l+d_rear)+dF1/2.+tpar[2];
  gMC->Gspos("YPF1", 2, "YMO1", 0., 0., dz, 0, "ONLY"); 

//
// pipe and heating jackets outside bellows
//
// left side
  cpar0[0]=(zvac1-abs_l+d_rear)/2;
  cpar0[1]=r_vacu+(abs_l-d_rear-abs_c)*TMath::Tan(theta_open1);
  cpar0[2]=r_abs +(abs_l-d_rear-abs_c)*TMath::Tan(theta_open1);
  cpar0[3]=cpar0[1]+2.*cpar0[0]*TMath::Tan(theta_open1);
  cpar0[4]=cpar0[2]+2.*cpar0[0]*TMath::Tan(theta_open1);
  gMC->Gsvolu("YV11", "CONE", idtmed[1758], cpar0, 5);
//
// insulation
  cpar[0]=cpar0[0];
  cpar[1]=cpar0[1]+d_tube;
  cpar[2]=cpar0[1]+d_tube+d_insu;
  cpar[3]=cpar0[3]+d_tube;
  cpar[4]=cpar0[3]+d_tube+d_insu;
  gMC->Gsvolu("YI11", "CONE", idtmed[1753], cpar, 5);
  gMC->Gspos("YI11", 1, "YV11", 0., 0., 0., 0, "ONLY"); 
//
// clearance
  cpar[1]=cpar0[2]-d_prot-d_free;
  cpar[2]=cpar0[2]-d_prot;
  cpar[3]=cpar0[4]-d_prot-d_free;
  cpar[4]=cpar0[4]-d_prot;
  gMC->Gsvolu("YP11", "CONE", idtmed[1755], cpar, 5);
  gMC->Gspos("YP11", 1, "YV11", 0., 0., 0., 0, "ONLY"); 
  
  dz=-dl+cpar0[0];
  gMC->Gspos("YV11", 1, "YMO1", 0., 0., dz, 0, "ONLY"); 

// right side
  cpar0[0]=(zvac4-zvac3)/2;
  cpar0[1]=r2-d_vacu;
  cpar0[2]=r2       ;
  cpar0[3]=cpar0[1]+2.*cpar0[0]*TMath::Tan(theta_open2);
  cpar0[4]=cpar0[2]+2.*cpar0[0]*TMath::Tan(theta_open2);
  gMC->Gsvolu("YV12", "CONE", idtmed[1758], cpar0, 5);
//
// insulation
  cpar[0]=cpar0[0];
  cpar[1]=cpar0[1]+d_tube;
  cpar[2]=cpar0[1]+d_tube+d_insu;
  cpar[3]=cpar0[3]+d_tube;
  cpar[4]=cpar0[3]+d_tube+d_insu;
  gMC->Gsvolu("YI12", "CONE", idtmed[1753], cpar, 5);
  gMC->Gspos("YI12", 1, "YV12", 0., 0., 0., 0, "ONLY"); 
//
// clearance
  cpar[1]=cpar0[2]-d_prot-d_free;
  cpar[2]=cpar0[2]-d_prot;
  cpar[3]=cpar0[4]-d_prot-d_free;
  cpar[4]=cpar0[4]-d_prot;
  gMC->Gsvolu("YP12", "CONE", idtmed[1755], cpar, 5);
  gMC->Gspos("YP12", 1, "YV12", 0., 0., 0., 0, "ONLY"); 
  
  dz=dl-cpar0[0];
  gMC->Gspos("YV12", 1, "YMO1", 0., 0., dz, 0, "ONLY"); 

//
// Second Section

  par2[0]  = 0.;
  par2[1]  = 360.;
  par2[2]  = 6.;
  dl=(zvac7-zvac4)/2.;
  
  par2[3]  = -dl;
  par2[4]  = r2+(zvac5-zvac3) * TMath::Tan(theta_open2);
  par2[5]  = R11;
  
  par2[6]  = -dl;
  par2[7]  = par2[4];
  par2[8]  = R21;

  par2[9]  = -dl+(zvac6-zvac4);
  par2[10]  = r2+(zvac6-zvac3) * TMath::Tan(theta_open2);
  par2[11]  = R21;

  par2[12] = -dl+(zvac6-zvac4);
  par2[13] = par2[10];
  par2[14] = zvac6*TMath::Tan(acc_min);

  par2[15] = -dl+(zcone_e-zvac4);
  par2[16] = r2+(zcone_e-zvac3) * TMath::Tan(theta_open2);
  par2[17] = 30.;

  par2[18] = -dl+(zvac7-zvac4);
  par2[19] = r2+(zvac7-zvac3) * TMath::Tan(theta_open2);
  par2[20] = 30.;

  
  gMC->Gsvolu("YGO2", "PCON", idtmed[1760], par2, 21);

  for (Int_t i=4; i<20; i+=3) par2[i]  = 0;
      
  gMC->Gsvolu("YMO2", "PCON", idtmed[1755], par2, 21);
  gMC->Gspos("YGO2", 1, "YMO2", 0., 0., 0., 0, "ONLY");  

  DZ+=dl;
  gMC->Gspos("YMO2", 1, "YMOT", 0., 0., DZ, 0, "ONLY");  
  DZ+=dl;
//
// Steel envelope
//
  tpar[0]=R11-dRSteel2;
  tpar[1]=R21;
  tpar[2]=2;
  gMC->Gsvolu("YS21", "TUBE", idtmed[1718], tpar, 3);
  dz=-dl+tpar[2];
  gMC->Gspos("YS21", 1, "YGO2", 0., 0., dz, 0, "ONLY");  
  dz+=tpar[2];
  tpar[0]=R21-dRSteel2;
  tpar[1]=R21;
  tpar[2]=(zvac6-zvac5)/2.;
  gMC->Gsvolu("YS22", "TUBE", idtmed[1718], tpar, 3);
  dz+=tpar[2];
  gMC->Gspos("YS22", 1, "YGO2", 0., 0., dz, 0, "ONLY");  
  dz+=tpar[2];
  
  cpar[0]=2.;
  cpar[1]=R21-dRSteel2;
  cpar[2]=zvac6 * TMath::Tan(acc_min);
  cpar[3]=cpar[1];
  cpar[4]=cpar[2]+4.*TMath::Tan(acc_min);

  gMC->Gsvolu("YS23", "CONE", idtmed[1718], cpar, 5);
  dz+=cpar[0];
  gMC->Gspos("YS23", 1, "YGO2", 0., 0., dz, 0, "ONLY");  
  dz+=cpar[0];

  cpar[0]=(zcone_e-zvac6-4.)/2;
  cpar[2]=cpar[4];
  cpar[4]=cpar[2]+2.*cpar[0]*TMath::Tan(acc_min);
  cpar[1]=cpar[2]-dRSteel2;
  cpar[3]=cpar[4]-dRSteel2;

  gMC->Gsvolu("YS24", "CONE", idtmed[1718], cpar, 5);
  dz+=cpar[0];
  gMC->Gspos("YS24", 1, "YGO2", 0., 0., dz, 0, "ONLY");  
  dz+=cpar[0];

  tpar[0]=26.;
  tpar[1]=30.;
  tpar[2]=(zvac7-zcone_e)/2.;

  gMC->Gsvolu("YS25", "TUBE", idtmed[1718], tpar, 3);
  dz+=tpar[2];
  gMC->Gspos("YS25", 1, "YGO2", 0., 0., dz, 0, "ONLY");  
  dz+=tpar[2];

//
// 2nd section: vacuum system 

  cpar0[0]=(zvac7-zvac4)/2;
  cpar0[1]=r2-d_vacu+(zvac4-zvac3)*TMath::Tan(theta_open2);
  cpar0[2]=r2       +(zvac4-zvac3)*TMath::Tan(theta_open2);
  cpar0[3]=cpar0[1]+2.*cpar0[0]*TMath::Tan(theta_open2);
  cpar0[4]=cpar0[2]+2.*cpar0[0]*TMath::Tan(theta_open2);
  gMC->Gsvolu("YV21", "CONE", idtmed[1758], cpar0, 5);
//
// insulation
  cpar[0]=cpar0[0];
  cpar[1]=cpar0[1]+d_tube;
  cpar[2]=cpar0[1]+d_tube+d_insu;
  cpar[3]=cpar0[3]+d_tube;
  cpar[4]=cpar0[3]+d_tube+d_insu;
  gMC->Gsvolu("YI21", "CONE", idtmed[1753], cpar, 5);
  gMC->Gspos("YI21", 1, "YV21", 0., 0., 0., 0, "ONLY"); 
//
// clearance
  cpar[1]=cpar0[2]-d_prot-d_free;
  cpar[2]=cpar0[2]-d_prot;
  cpar[3]=cpar0[4]-d_prot-d_free;
  cpar[4]=cpar0[4]-d_prot;
  gMC->Gsvolu("YP21", "CONE", idtmed[1755], cpar, 5);
  gMC->Gspos("YP21", 1, "YV21", 0., 0., 0., 0, "ONLY"); 
  
  dz=0.;
  gMC->Gspos("YV21", 1, "YMO2", 0., 0., dz, 0, "ONLY"); 

//
// Third Section: Bellows and Flange 
//
  par3[0]  = 0.;
  par3[1]  = 360.;
  par3[2]  = 8.;
  dl=(zvac9-zvac7)/2.;
  
  par3[3]  = -dl;
  par3[4]  = r2+(zvac7-zvac3) * TMath::Tan(theta_open2);
  par3[5]  = 30.;

  par3[6]  = -dl+dr21;
  par3[7]  = par3[4]+dr21;
  par3[8]  = 30.;

  par3[9]  = par3[6]+dB2;
  par3[10] = par3[7];
  par3[11] = 30.;

  par3[12] = par3[9]+dr22;
  par3[13] = par3[10]+dr22;
  par3[14] = 30.;

  par3[15] = par3[12]+dF2;
  par3[16] = par3[13];
  par3[17] = 30.;

  par3[18] = par3[15]+dr22;
  par3[19] = par3[16]-dr22;
  par3[20] = 30.;

  par3[21] = par3[18]+dB2;
  par3[22] = par3[19];
  par3[23] = 30.;

  par3[24] = par3[21]+dr23;
  par3[25] = par3[22]-dr23;
  par3[26] = 30.;

  Float_t r3=par3[25];
  
  gMC->Gsvolu("YGO3", "PCON", idtmed[1760], par3, 27);

  for (Int_t i=4; i<26; i+=3) par3[i]  = 0;
      
  gMC->Gsvolu("YMO3", "PCON", idtmed[1755], par3, 27);
  gMC->Gspos("YGO3", 1, "YMO3", 0., 0., 0., 0, "ONLY");  

//
// Steel envelope
  tpar[0]=26;
  tpar[1]=30;
  tpar[2]=dl;
  gMC->Gsvolu("YS31", "TUBE", idtmed[1718], tpar, 3);
  gMC->Gspos("YS31", 1, "YGO3", 0., 0., 0., 0, "ONLY");  
  DZ+=dl;
  gMC->Gspos("YMO3", 1, "YMOT", 0., 0., DZ, 0, "ONLY");  
  DZ+=dl;

//
// 3rd section: vacuum system
//
//
// Bellow2
//
  tpar[0]=rB2;
  tpar[1]=rB2+hB2;
  tpar[2]=eB2/2.;
  gMC->Gsvolu("YB21", "TUBE", idtmed[1758], tpar, 3);
  dl1=tpar[2];
  
  tpar[0]=rB2+hB2-eB2;
  tpar[1]=rB2+hB2;
  tpar[2]=(lB2/2.-2.*eB2)/2.;
  gMC->Gsvolu("YB22", "TUBE", idtmed[1758], tpar, 3);
  dl2=tpar[2];

  tpar[0]=rB2-eB2;
  tpar[1]=rB2;
  tpar[2]=lB2/8.;
  gMC->Gsvolu("YB23", "TUBE", idtmed[1758], tpar, 3);
  dl3=tpar[2];


  tpar[0]=0;
  tpar[1]=rB2+hB2;
  tpar[2]=lB2/2.;
  gMC->Gsvolu("YBU2", "TUBE", idtmed[1755], tpar, 3);

  dz=-tpar[2]+dl3;
  gMC->Gspos("YB23", 1, "YBU2", 0., 0., dz, 0, "ONLY"); 
  dz+=dl3;
  dz+=dl1;  
  gMC->Gspos("YB21", 1, "YBU2", 0., 0., dz, 0, "ONLY"); 
  dz+=dl1;  
  dz+=dl2;  
  gMC->Gspos("YB22", 1, "YBU2", 0., 0., dz, 0, "ONLY"); 
  dz+=dl2;  
  dz+=dl1;
  gMC->Gspos("YB21", 2, "YBU2", 0., 0., dz, 0, "ONLY"); 
  dz+=dl1;
  dz+=dl3;
  gMC->Gspos("YB23", 2, "YBU2", 0., 0., dz, 0, "ONLY"); 
  

  tpar[0]=0;
  tpar[1]=rB2+hB2;
  tpar[2]=7.*lB2/2.;
  gMC->Gsvolu("YBM2", "TUBE", idtmed[1755], tpar, 3);
  dz=-tpar[2]+lB2/2.;
  for (Int_t i=0; i<7; i++) {
      gMC->Gspos("YBU2", i+1 , "YBM2", 0., 0.,dz , 0, "ONLY"); 
      dz+=lB2;
  }

  dz=-dl+dr21+tpar[2];
  gMC->Gspos("YBM2", 1, "YMO3", 0., 0., dz, 0, "ONLY"); 

  dz=dl-dr23-tpar[2];
  gMC->Gspos("YBM2", 2, "YMO3", 0., 0., dz, 0, "ONLY"); 

//
// Flange

  tpar[0]=0;
  tpar[1]=rF2;
  tpar[2]=dF2/2.;
  gMC->Gsvolu("YFM2", "TUBE", idtmed[1755], tpar, 3);

  tpar[0]=rF2-d_flange;
  tpar[1]=rF2;
  tpar[2]=dF2/2.;
  gMC->Gsvolu("YF21", "TUBE", idtmed[1758], tpar, 3);
  gMC->Gspos("YF21", 1, "YFM2", 0., 0., 0., 0, "ONLY"); 

  tpar[0]=rB2;
  tpar[1]=rF2-d_flange;
  tpar[2]=d_flange/2.;
  gMC->Gsvolu("YF22", "TUBE", idtmed[1758], tpar, 3);
  dz=-dF2/2.+tpar[2];
  gMC->Gspos("YF22", 1, "YFM2", 0., 0., dz, 0, "ONLY"); 
  dz= dF2/2.-tpar[2];
  gMC->Gspos("YF22", 2, "YFM2", 0., 0., dz, 0, "ONLY"); 

  dz=dr21/2.-dr23/2.;
  gMC->Gspos("YFM2", 2, "YMO3", 0., 0., dz, 0, "ONLY"); 


//
// pipe between flange and bellows
  tpar[0]=rB2-d_tube;
  tpar[1]=rB2;
  tpar[2]=2.*(dB2+dr22-7.*lB2)/4.;
  gMC->Gsvolu("YPF2", "TUBE", idtmed[1758], tpar, 3);
  dz=dr21/2.-dr23/2.-dF2/2.-tpar[2];
  gMC->Gspos("YPF2", 1, "YMO3", 0., 0., dz, 0, "ONLY"); 
  dz=dr21/2.-dr23/2.+dF2/2.+tpar[2];
  gMC->Gspos("YPF2", 2, "YMO3", 0., 0., dz, 0, "ONLY"); 

//
// 4th section: rear shield and closing cone
//
  par4[0]  = 0.;
  par4[1]  = 360.;
  par4[2]  = 6.;
  dl=(zvac12-zvac9)/2.;
  
  par4[3]  = -dl;
  par4[4]  = r3;
  par4[5]  = 30.;

  par4[6]  = -dl+(zvac10-zvac9);
  par4[7]  = r3+(zvac10-zvac9) * TMath::Tan(theta_open3);
  par4[8]  = 30.;

  par4[9]  = par4[6];
  par4[10] = par4[7];
  par4[11] = R42;

  par4[12] = -dl+(zvac11-zvac9);
  par4[13] = r3+(zvac11-zvac9) * TMath::Tan(theta_open3);
  par4[14] = R42;

  par4[15] = par4[12];
  par4[16] = par4[13];
  par4[17] = R43;

  par4[18] = -dl+(zvac12-zvac9);
  par4[19] = r_abs;
  par4[20] = R43;

  gMC->Gsvolu("YGO4", "PCON", idtmed[1760], par4, 21);
  for (Int_t i=4; i<20; i+=3) par4[i]  = 0;
      

  gMC->Gsvolu("YMO4", "PCON", idtmed[1755], par4, 21);
  gMC->Gspos("YGO4", 1, "YMO4", 0., 0., 0., 0, "ONLY");  


  DZ+=dl;
  gMC->Gspos("YMO4", 1, "YMOT", 0., 0., DZ, 0, "ONLY");  
  DZ+=dl;
//
// Closing concrete cone 
//
  cpar[0]=(zvac12-zvac11)/2.;
  cpar[1] = r3+(zvac11-zvac9) * TMath::Tan(theta_open3);
  cpar[2] = cpar[1]+0.001;
  cpar[3] = r_abs;
  cpar[4] = cpar[2];
  gMC->Gsvolu("YCC4", "CONE", idtmed[1752], cpar, 5);
  dz=dl-cpar[0];
  gMC->Gspos("YCC4", 1, "YGO4", 0., 0., dz, 0, "ONLY");  
//
// Steel envelope
//
  dz=-dl;
  tpar[0]=26.;
  tpar[1]=30.;
  tpar[2]=(zvac10-zvac9)/2.;
  gMC->Gsvolu("YS41", "TUBE", idtmed[1718], tpar, 3);
  dz+=tpar[2];
  gMC->Gspos("YS41", 1, "YGO4", 0., 0., dz, 0, "ONLY");  
  dz+=tpar[2];

  tpar[0]=26.;
  tpar[1]=R41;
  tpar[2]=2.;
  gMC->Gsvolu("YS42", "TUBE", idtmed[1718], tpar, 3);
  dz+=tpar[2];
  gMC->Gspos("YS42", 1, "YGO4", 0., 0., dz, 0, "ONLY");  
  dz+=tpar[2];

  tpar[0]=R41-dRSteel2;
  tpar[1]=R41;
  tpar[2]=(zvac11-zvac10-4)/2.;
  gMC->Gsvolu("YS43", "TUBE", idtmed[1718], tpar, 3);
  dz+=tpar[2];
  gMC->Gspos("YS43", 1, "YGO4", 0., 0., dz, 0, "ONLY");  
//
// rear lead shield
//
  tpar[0]=R41;
  tpar[1]=R42;
  tpar[2]=(zvac11-zvac10)/2.;
  gMC->Gsvolu("YPBI", "TUBE", idtmed[1752], tpar, 3);
  dz-=4;
  gMC->Gspos("YPBI", 1, "YGO4", 0., 0., dz, 0, "ONLY"); 

  tpar[0]=R42-5;
  tpar[1]=R42;
  tpar[2]=(zvac11-zvac10)/2.;
  gMC->Gsvolu("YPBO", "TUBE", idtmed[1712], tpar, 3);
  gMC->Gspos("YPBO", 1, "YPBI", 0., 0., 0., 0, "ONLY"); 
  
//
// rear Fe shield
//

  tpar[0]=R42;
  tpar[1]=R43;
  tpar[2]=(zvac12-zvac11)/2.;
  gMC->Gsvolu("YFEI", "TUBE", idtmed[1749], tpar, 3);
  dz=dl-tpar[2];
  gMC->Gspos("YFEI", 1, "YGO4", 0., 0., dz, 0, "ONLY"); 

  tpar[0]=R42;
  tpar[1]=R43;
  tpar[2]=2.5;
  gMC->Gsvolu("YFEO", "TUBE", idtmed[1709], tpar, 3);
  dz=-(zvac12-zvac11)/2.+tpar[2];
  gMC->Gspos("YFEO", 1, "YFEI", 0., 0., dz, 0, "ONLY"); 

//
// 4th section: vacuum system 
//
// up to closing cone
  cpar0[0]=(zvac11-zvac9)/2;
  cpar0[1]=r3-d_vacu;
  cpar0[2]=r3;
  cpar0[3]=cpar0[1]+2.*cpar0[0]*TMath::Tan(theta_open3);
  cpar0[4]=cpar0[2]+2.*cpar0[0]*TMath::Tan(theta_open3);
  gMC->Gsvolu("YV31", "CONE", idtmed[1758], cpar0, 5);
//
// insulation
  cpar[0]=cpar0[0];
  cpar[1]=cpar0[1]+d_tube;
  cpar[2]=cpar0[1]+d_tube+d_insu;
  cpar[3]=cpar0[3]+d_tube;
  cpar[4]=cpar0[3]+d_tube+d_insu;
  gMC->Gsvolu("YI31", "CONE", idtmed[1753], cpar, 5);
  gMC->Gspos("YI31", 1, "YV31", 0., 0., 0., 0, "ONLY"); 
//
// clearance
  cpar[1]=cpar0[2]-d_prot-d_free;
  cpar[2]=cpar0[2]-d_prot;
  cpar[3]=cpar0[4]-d_prot-d_free;
  cpar[4]=cpar0[4]-d_prot;
  gMC->Gsvolu("YP31", "CONE", idtmed[1755], cpar, 5);
  gMC->Gspos("YP31", 1, "YV31", 0., 0., 0., 0, "ONLY"); 
  
  dz=-dl+cpar[0];
  gMC->Gspos("YV31", 1, "YMO4", 0., 0., dz, 0, "ONLY"); 
//
// closing cone
  cpar0[0]=(zvac12-zvac11)/2;
  cpar0[1]=r3-d_vacu+(zvac11-zvac9)*TMath::Tan(theta_open3);
  cpar0[2]=r3       +(zvac11-zvac9)*TMath::Tan(theta_open3);
  cpar0[3]=r_vacu;
  cpar0[4]=r_abs;
  gMC->Gsvolu("YV32", "CONE", idtmed[1758], cpar0, 5);
//
// insulation
  cpar[0]=cpar0[0];
  cpar[1]=cpar0[1]+d_tube;
  cpar[2]=cpar0[1]+d_tube+d_insu;
  cpar[3]=cpar0[3]+d_tube;
  cpar[4]=cpar0[3]+d_tube+d_insu;
  gMC->Gsvolu("YI32", "CONE", idtmed[1753], cpar, 5);
  gMC->Gspos("YI32", 1, "YV32", 0., 0., 0., 0, "ONLY"); 
//
// clearance
  cpar[1]=cpar0[2]-d_prot-d_free;
  cpar[2]=cpar0[2]-d_prot;
  cpar[3]=cpar0[4]-d_prot-d_free;
  cpar[4]=cpar0[4]-d_prot;
  gMC->Gsvolu("YP32", "CONE", idtmed[1755], cpar, 5);
  gMC->Gspos("YP32", 1, "YV32", 0., 0., 0., 0, "ONLY"); 
  
  dz=dl-cpar[0];
  gMC->Gspos("YV32", 1, "YMO4", 0., 0., dz, 0, "ONLY"); 
//
//
// MUON trigger wall
// 
  tpar[0] = 30.;
  tpar[1] = 310.;
  tpar[2] = (zfil_out - zfil_in) / 2.;
  gMC->Gsvolu("YFIM", "TUBE", idtmed[1749], tpar, 3);
  dz = (zfil_in + zfil_out) / 2.;
  tpar[2] -= 10.;
  gMC->Gsvolu("YFII","TUBE", idtmed[1709], tpar, 3);
  gMC->Gspos("YFII", 1, "YFIM", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("YFIM", 1, "ALIC", 0., 0., dz, 0, "ONLY");
//
// Shielding close to chamber
//
  
  cpar[0]=(zch1-dzch-1.-abs_l)/2.;
  cpar[1]=R11;
  cpar[2]=abs_l*TMath::Tan(acc_min);
  cpar[3]=R11;
  cpar[4]=(abs_l+2.*cpar[0])*TMath::Tan(acc_min);
  gMC->Gsvolu("YCS1", "CONE", idtmed[1720], cpar, 5);
  dz=abs_l+cpar[0];
  gMC->Gspos("YCS1", 1, "ALIC", 0., 0., dz, 0, "ONLY");

  cpar[0]=(zvac4-(zch1+dzch+1.))/2.;
  cpar[1]=R11;
  cpar[2]=(zvac4-2.*cpar[0])*TMath::Tan(acc_min);
  cpar[3]=R11;
  cpar[4]=zvac4*TMath::Tan(acc_min);
  gMC->Gsvolu("YCS2", "CONE", idtmed[1720], cpar, 5);
  dz=zvac4-cpar[0];
  gMC->Gspos("YCS2", 1, "ALIC", 0., 0., dz, 0, "ONLY");


  cpar[0]=(zch2-dzch-1.-zvac4)/2.;
  cpar[1]=R21;
  cpar[2]=zvac4*TMath::Tan(acc_min);
  cpar[3]=R21;
  cpar[4]=(zvac4+2.*cpar[0])*TMath::Tan(acc_min);
  gMC->Gsvolu("YCS3", "CONE", idtmed[1720], cpar, 5);
  dz=zvac4+cpar[0];
  gMC->Gspos("YCS3", 1, "ALIC", 0., 0., dz, 0, "ONLY");
  

  cpar[0]=(zvac6-(zch2+dzch+1.))/2.;
  cpar[1]=R11;
  cpar[2]=(zvac6-2.*cpar[0])*TMath::Tan(acc_min);
  cpar[3]=R11;
  cpar[4]=zvac6*TMath::Tan(acc_min);
  gMC->Gsvolu("YCS4", "CONE", idtmed[1720], cpar, 5);
  dz=zvac6-cpar[0];
  gMC->Gspos("YCS4", 1, "ALIC", 0., 0., dz, 0, "ONLY");

}

void AliSHILv0::Init()
{
  //
  // Initialise the muon shield after it has been built
  //
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" SHILv0_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  // Here the SHIL initialisation code (if any!)
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

 







