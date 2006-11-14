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

#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliRun.h"
#include "AliSHILv0.h"
#include "AliLog.h"

ClassImp(AliSHILv0)
 
//_____________________________________________________________________________
AliSHILv0::AliSHILv0():
    fPbCone(1)
{
  //
  // Default constructor for muon shield
  //
    
}
 
//_____________________________________________________________________________
AliSHILv0::AliSHILv0(const char *name, const char *title)
    : AliSHIL(name,title), 
      fPbCone(1)
{
  //
  // Standard constructor for muon shield
  //
  // Pb  cone not yet compatible with muon chamber inner radii
  // Switched off by default
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

    Float_t cpar[5], cpar0[5], tpar[3], par1[39], pars1[100], par2[36], par3[27], 
	par4[24], par0[87];
    Float_t dz, dZ;
  
    Int_t *idtmed = fIdtmed->GetArray()-1699;

#include "ABSOSHILConst.h"
#include "SHILConst.h"

enum {kC=1705, kAl=1708, kFe=1709, kCu=1710, kW=1711, kPb=1712,
		kNiCuW=1720, kVacuum=1715, kAir=1714, kConcrete=1716,
		kPolyCH2=1717, kSteel=1709, kInsulation=1713};	
//
// Material of the rear part of the shield
  Int_t iHeavy=kNiCuW;
  if (fPbCone) iHeavy=kPb;
//
// Mother volume
//
  Float_t dRear1=kDRear;
  
  Float_t zstart=kZRear-dRear1;
  
  par0[0]  = 0.;
  par0[1]  = 360.;
  par0[2]  = 28.;

  Float_t dl=(kZvac12-zstart)/2.;
  dz=zstart+dl;
//
// start
  par0[3]  = -dl;
  par0[4]  = 0.;
  par0[5]  = zstart * TMath::Tan(kAccMin);
// recess station 1
  par0[6]  = -dz+kZch11;
  par0[7]  = 0.;
  par0[8]  = kZch11 * TMath::Tan(kAccMin);

  par0[9]   = par0[6];
  par0[10]  = 0.;
  par0[11]  = 17.9;

  par0[12]  = -dz+kZch12;
  par0[13]  = 0.;
  par0[14]  = 17.9;

  par0[15]  = par0[12];
  par0[16]  = 0.;
  par0[17]  = kZch12 * TMath::Tan(kAccMin);
// recess station 2
  par0[18]  = -dz+kZch21;
  par0[19]  = 0.;
  par0[20]  = kZch21 * TMath::Tan(kAccMin);

  par0[21]  = -dz+kZch21;
  par0[22] = 0.;
  par0[23] = 23.;

  par0[24]  = -dz+kZch22;
  par0[25] = 0.;
  par0[26] = 23.;

  par0[27]  = -dz+kZch22;
  par0[28]  = 0.;
  par0[29]  = kZch22 * TMath::Tan(kAccMin);
//
  par0[30] = -dz+kZvac6;
  par0[31] = 0.;
  par0[32] = kZvac6 * TMath::Tan(kAccMin);
// end of 2 deg cone
  par0[33] = -dz+kZConeE;
  par0[34] = 0.;
  par0[35] = 30.;

  par0[36] = -dz+kZch31;
  par0[37] = 0.;
  par0[38] = 30.;

  par0[39] = -dz+kZch31;
  par0[40] = 0.;
  par0[41] = 29.;

  par0[42] = -dz+kZch32;
  par0[43] = 0.;
  par0[44] = 29.;
// start of 1.6 deg cone
  par0[45] = -dz+kZch32;
  par0[46] = 0.;
  par0[47] = 30.+(kZch32-kZConeE)*TMath::Tan(kThetaOpenPbO);
// recess station 4
  par0[48] = -dz+kZch41;
  par0[49] = 0.;
  par0[50] = 30.+(kZch41-kZConeE)*TMath::Tan(kThetaOpenPbO);

  par0[51] = -dz+kZch41;
  par0[52] = 0.;
  par0[53] = 37.5;

  par0[54] = -dz+kZch42;
  par0[55] = 0.;
  par0[56] = 37.5;

  par0[57] = -dz+kZch42;
  par0[58] = 0.;
  par0[59] = 30.+(kZch42-kZConeE)*TMath::Tan(kThetaOpenPbO);

// recess station 5

  par0[60] = -dz+kZch51;
  par0[61] = 0.;
  par0[62] = 30.+(kZch51-kZConeE)*TMath::Tan(kThetaOpenPbO);

  par0[63] = -dz+kZch51;
  par0[64] = 0.;
  par0[65] = 37.5;

  par0[66] = -dz+kZch52;
  par0[67] = 0.;
  par0[68] = 37.5;

  par0[69] = -dz+kZch52;
  par0[70] = 0.;
  par0[71] = 30.+(kZch52-kZConeE)*TMath::Tan(kThetaOpenPbO);

// end of cone

  par0[72] = -dz+kZvac10;
  par0[73] = 0.;
  par0[74] = 30.+(kZvac10-kZConeE)*TMath::Tan(kThetaOpenPbO);

  par0[75] = -dz+kZvac10;
  par0[76] = 0.;
  par0[77] = kR42;

  par0[78] = -dz+kZvac11;
  par0[79] = 0.;
  par0[80] = kR42;

  par0[81] = -dz+kZvac11;
  par0[82] = 0.;
  par0[83] = kR43;

  par0[84] = -dz+kZvac12;
  par0[85] = 0.;
  par0[86] = kR43;

  gMC->Gsvolu("YMOT", "PCON", idtmed[kVacuum], par0, 87);
  dz=zstart+dl;
  gMC->Gspos("YMOT", 1, "ALIC", 0., 0., dz, 0, "ONLY");  
  gMC->Gsbool("YMOT","L3DO");
  gMC->Gsbool("YMOT","L3O1");
  gMC->Gsbool("YMOT","L3O2");
//

  dZ=-dl;

//
// First section: bellows below and behind front absorber 
// 
//
  par1[0]  = 0.;
  par1[1]  = 360.;
  par1[2]  = 12.;
  dl=(kZvac4-zstart)/2.;
  
  par1[3]  = -dl;
  par1[4]  = kRAbs+(zstart-kZOpen) * TMath::Tan(kThetaOpen1);
  par1[5]  = zstart * TMath::Tan(kAccMin);

  par1[6]  = -dl+kZvac1-zstart;
  par1[7]  = kRAbs+ (kZvac1-kZOpen) * TMath::Tan(kThetaOpen1);
  par1[8]  = kZvac1 * TMath::Tan(kAccMin);

  par1[9]  = par1[6]+kDr11/2.;
  par1[10] = par1[7]+kDr11;
  par1[11] = (kZvac1+kDr11/2.) * TMath::Tan(kAccMin);

  par1[12] = -dl+dRear1;
  par1[13] = par1[10];
  par1[14] = kZRear * TMath::Tan(kAccMin);

  par1[15] = -dl+dRear1;
  par1[16] = par1[10];
  par1[17] = kR11;

  par1[18] = -dl+(kZvac1+kDr11+kDB1-zstart);
  par1[19] = par1[16];
  par1[20] = kR11;

  par1[21] = par1[18]+kDr12;
  par1[22] = par1[19]+kDr12;
  par1[23] = kR11;

  par1[24] = par1[21]+kDF1;
  par1[25] = par1[22];
  par1[26] = kR11;

  par1[27] = par1[24]+kDr12;
  par1[28] = par1[25]-kDr12; 
  par1[29] = kR11;

  par1[30] = par1[27]+kDB1;
  par1[31] = par1[28];
  par1[32] = kR11;

  par1[33] = par1[30]+kDr13;
  par1[34] = par1[31]-kDr13;
  par1[35] = kR11;

  par1[36] = -dl+kZvac4-zstart;
  par1[37] = par1[34];
  par1[38] = kR11;

  Float_t r2  = par1[37];
  Float_t rBox= par1[31]-0.1;

  gMC->Gsvolu("YGO1", "PCON", idtmed[kNiCuW], par1, 39);
  Int_t i;
  
  for (i=0; i<39; i++)  pars1[i]  = par1[i];
  for (i=4; i<38; i+=3) pars1[i]  = 0.;

  gMC->Gsvolu("YMO1", "PCON", idtmed[kVacuum+40], pars1, 39);
  gMC->Gspos("YGO1", 1, "YMO1", 0., 0., 0., 0, "ONLY");  
  dZ+=dl;
  gMC->Gspos("YMO1", 1, "YMOT", 0., 0., dZ, 0, "ONLY");  
  dZ+=dl;

//
// Steel envelope
  tpar[0]=kR11-kDRSteel2;
  tpar[1]=kR11;
  tpar[2]=(kZvac4-kZvac3)/2.;
  gMC->Gsvolu("YSE1", "TUBE", idtmed[kNiCuW], tpar, 3);
  dz=dl-tpar[2];
  gMC->Gspos("YSE1", 1, "YGO1", 0., 0., dz, 0, "ONLY");

//
// 1st section: vacuum system
//
//
// Bellow 1
//

//
// Bellow 1
//
  tpar[0]=kRB1;
  tpar[1]=kRB1+kHB1;
  tpar[2]=kEB1/2.;
  gMC->Gsvolu("YB11", "TUBE", idtmed[kSteel+40], tpar, 3);
  Float_t dl1=tpar[2];
  
  tpar[0]=kRB1+kHB1-kEB1;
  tpar[1]=kRB1+kHB1;
  tpar[2]=(kLB1/2.-2.*kEB1)/2.;
  gMC->Gsvolu("YB12", "TUBE", idtmed[kSteel+40], tpar, 3);
  Float_t dl2=tpar[2];

  tpar[0]=kRB1-kEB1;
  tpar[1]=kRB1;
  tpar[2]=kLB1/8.;
  gMC->Gsvolu("YB13", "TUBE", idtmed[kSteel+40], tpar, 3);
  Float_t dl3=tpar[2];


  tpar[0]=0;
  tpar[1]=kRB1+kHB1;
  tpar[2]=-kLB1/2.;
  gMC->Gsvolu("YBU1", "TUBE", idtmed[kVacuum+40], tpar, 3);

  dz=-kLB1/2.+dl3;
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
  tpar[1]=kRB1+kHB1+0.5;
  tpar[2]=12.*kLB1/2.;
  gMC->Gsvolu("YBM1", "TUBE", idtmed[kVacuum+40], tpar, 3);
  gMC->Gsdvn("YB1S", "YBM1", 12 , 3);

  Float_t bsize = tpar[2];
  tpar[0]=kRB1+kHB1;
  tpar[2]=-kLB1/2.;
  gMC->Gsvolu("YBI1", "TUBE", idtmed[kInsulation+40], tpar, 3);

  gMC->Gspos("YBI1", 1, "YB1S", 0., 0., 0., 0, "ONLY"); 
  gMC->Gspos("YBU1", 1, "YB1S", 0., 0., 0., 0, "ONLY"); 

  dz=-dl+(kZvac1-zstart)+kDr11/2.+bsize;
  gMC->Gspos("YBM1", 1, "YMO1", 0., 0., dz, 0, "ONLY"); 

//
// Flange

  tpar[0]=0;
  tpar[1]=kRF1+0.6;
  tpar[2]=kDF1/2.;
  gMC->Gsvolu("YFM1", "TUBE", idtmed[kVacuum+40], tpar, 3);
// Steel
  tpar[0]=kRB1;
  tpar[1]=kRF1+0.6;
  tpar[2]=kDF1/2.;
  gMC->Gsvolu("YF11", "TUBE", idtmed[kSteel+40], tpar, 3);
// Insulation
  tpar[0]=kRF1;
  tpar[1]=kRF1+0.5;
  tpar[2]=kDF1/2.;
  gMC->Gsvolu("YF12", "TUBE", idtmed[kInsulation+40], tpar, 3);


  gMC->Gspos("YF11", 1, "YFM1", 0., 0., 0., 0, "ONLY"); 
  gMC->Gspos("YF12", 1, "YFM1", 0., 0., 0., 0, "ONLY"); 

  dz=-dl+(kZvac1-zstart)+kDr11/2.+2.*bsize+kDF1/2.+3.;
  gMC->Gspos("YFM1", 2, "YMO1", 0., 0., dz, 0, "ONLY"); 

//
// pipe between flange and bellows
//
// Steel 
  tpar[0]=kRB1-dTubeS;
  tpar[1]=kRB1+0.6;
  tpar[2]=1.5;
  gMC->Gsvolu("YPF1", "TUBE", idtmed[kSteel+40], tpar, 3);
// Insulation
  tpar[0]=kRB1;
  tpar[1]=kRB1+0.5;
  gMC->Gsvolu("YPS1", "TUBE", idtmed[kInsulation+40], tpar, 3);
  gMC->Gspos("YPS1", 1, "YPF1", 0., 0., 0., 0, "ONLY"); 

  dz=dz-1.5-kDF1/2.;
  gMC->Gspos("YPF1", 1, "YMO1", 0., 0., dz, 0, "ONLY"); 
  dz=dz+3.0+kDF1;
  gMC->Gspos("YPF1", 2, "YMO1", 0., 0., dz, 0, "ONLY"); 
//

// Pipe+Heating     1.5 mm 
// Heating Jacket   5.0 mm
// Protection       1.0 mm
// ========================
//                  7.5 mm
// pipe and heating jackets outside bellows
//
// left side
  cpar0[0]=(kZvac1+kDr11/2.-zstart)/2;
  cpar0[1]=kRVacu-0.05  +(zstart-kZOpen)*TMath::Tan(kThetaOpen1);
  cpar0[2]=kRVacu+0.7   +(zstart-kZOpen)*TMath::Tan(kThetaOpen1);
  cpar0[3]=cpar0[1]+2.*cpar0[0]*TMath::Tan(kThetaOpen1);
  cpar0[4]=cpar0[2]+2.*cpar0[0]*TMath::Tan(kThetaOpen1);
  gMC->Gsvolu("YV11", "CONE", idtmed[kSteel+40], cpar0, 5);
//
// insulation
  dTubeS=0.15;
  cpar[0]=cpar0[0];
  cpar[1]=cpar0[1]+0.15;
  cpar[2]=cpar0[1]+0.65;
  cpar[3]=cpar0[3]+0.15;
  cpar[4]=cpar0[3]+0.65;
  gMC->Gsvolu("YI11", "CONE", idtmed[kInsulation+40], cpar, 5);
  gMC->Gspos("YI11", 1, "YV11", 0., 0., 0., 0, "ONLY"); 
  dz=-dl+cpar0[0];
  gMC->Gspos("YV11", 1, "YMO1", 0., 0., dz, 0, "ONLY"); 

// right side
  dTubeS  = 0.35;
  dVacuS += 0.25;
  
  cpar0[0] = (kZvac4-kZvac3)/2;
  cpar0[1] = kRB1;
  cpar0[2] = cpar0[1]+dVacuS;
  cpar0[3] = cpar0[1]+2.*cpar0[0]*TMath::Tan(kThetaOpenB);
  cpar0[4] = cpar0[2]+2.*cpar0[0]*TMath::Tan(kThetaOpenB);
  gMC->Gsvolu("YV12", "CONE", idtmed[kSteel], cpar0, 5);
  Float_t r2V=cpar0[3];
//
// insulation
  cpar[0] = cpar0[0];
  cpar[1] = cpar0[1]+dTubeS;
  cpar[2] = cpar0[1]+dTubeS+kDInsuS;
  cpar[3] = cpar0[3]+dTubeS;
  cpar[4] = cpar0[3]+dTubeS+kDInsuS;
  gMC->Gsvolu("YI12", "CONE", idtmed[kInsulation], cpar, 5);
  gMC->Gspos("YI12", 1, "YV12", 0., 0., 0., 0, "ONLY"); 

  dz=dl-cpar0[0];
  gMC->Gspos("YV12", 1, "YMO1", 0., 0., dz, 0, "ONLY"); 

//
// Second Section
// Between first and second bellow section
//

  par2[0]  = 0.;
  par2[1]  = 360.;
  par2[2]  = 11.;
  dl=(kZvac7-kZvac4)/2.;
// recess station 2
  par2[3]  = -dl;
  par2[4]  = r2;
  par2[5]  = kR21;

  par2[6]  = -dl+.1;
  par2[7]  = r2;
  par2[8]  = kR21;

  par2[9]   = -dl+(kZvac6-kZvac4);
  par2[10]  = r2+(kZvac6-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  par2[11]  = kR21;

  par2[12] = -dl+(kZvac6-kZvac4);
  par2[13] = par2[10];
  par2[14] = kZvac6*TMath::Tan(kAccMin);

// Start of Pb section
  par2[15] = -dl+(kZPb-kZvac4);
  par2[16] = r2+(kZPb-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  par2[17] = kZPb*TMath::Tan(kAccMin);
//
// end of cone following 2 deg line
  par2[18] = -dl+(kZConeE-kZvac4);
  par2[19] = r2+(kZConeE-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  par2[20] = 30.;
// recess station 3
  par2[21] = -dl+(kZch31-kZvac4);
  par2[22] = r2+(kZch31-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  par2[23] = 30.;

  par2[24] = -dl+(kZch31-kZvac4);
  par2[25] = r2+(kZch31-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  par2[26] = 29.;

  par2[27] = -dl+(kZch32-kZvac4);
  par2[28] = r2+(kZch32-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  par2[29] = 29.;

  par2[30] = -dl+(kZch32-kZvac4);
  par2[31] = r2+(kZch32-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  par2[32] = 30.;

  par2[33] = -dl+(kZvac7-kZvac4);
  par2[34] = r2+(kZvac7-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  par2[35] = 30.;

  gMC->Gsvolu("YGO2", "PCON", idtmed[kSteel+40], par2, 36);

//
// Lead cone 
//
  Float_t parPb[12];
  parPb[0]  = 0.;
  parPb[1]  = 360.;
  parPb[2]  = 3.;
  Float_t dlPb=(kZvac7-kZPb)/2.;
  
  parPb[3]  = -dlPb;
  parPb[4]  =  r2+(kZPb-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  parPb[5]  =  kZPb*TMath::Tan(kAccMin)-kDRSteel2;
  
  parPb[6]  = -dlPb+(kZConeE-kZPb);
  parPb[7]  =  r2+(kZConeE-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  parPb[8]  = 26.;
  
  parPb[9]   = dlPb;
  parPb[10]  =  r2+(kZvac7-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  parPb[11]  = 26.;

  gMC->Gsvolu("YXO2", "PCON", idtmed[kPb], parPb, 12);	  
  gMC->Gspos("YXO2", 1, "YGO2", 0., 0., (kZPb-kZvac4)/2., 0, "ONLY");  

//
// W cone 
//
  Float_t parW[15];
  parW[0]  = 0.;
  parW[1]  = 360.;
  parW[2]  = 4.;
  Float_t dlW=(kZPb-kZvac4)/2.;
  
  parW[3]   = -dlW;
  parW[4]   =  r2;
  parW[5]   =  kR21-kDRSteel2;
  
  parW[6]   = -dlW+(kZvac6-kZvac4)+kDRSteel2;
  parW[7]   =  r2+(kZvac6-kZvac4+kDRSteel2) * TMath::Tan(kThetaOpen2);
  parW[8]   =  kR21-kDRSteel2;
 
  parW[9]   = -dlW+(kZvac6-kZvac4)+kDRSteel2;
  parW[10]  =  r2+(kZvac6-kZvac4+kDRSteel2) * TMath::Tan(kThetaOpen2);
  parW[11]  =  (kZvac6+kDRSteel2)*TMath::Tan(kAccMin)-kDRSteel2;
 
  parW[12]  = dlW;
  parW[13]  =  r2+(kZPb-kZvac4) * TMath::Tan(kThetaOpen2);
  parW[14]  = kZPb*TMath::Tan(kAccMin)-kDRSteel2;

  gMC->Gsvolu("YYO2", "PCON", idtmed[kNiCuW], parW, 15);	  
  gMC->Gspos("YYO2", 1, "YGO2", 0., 0., -(kZvac7-kZPb)/2., 0, "ONLY");  

  for (i=4; i<35; i+=3) par2[i]  = 0;
          
  gMC->Gsvolu("YMO2", "PCON", idtmed[kVacuum+40], par2, 36);
  gMC->Gspos("YGO2", 1, "YMO2", 0., 0., 0., 0, "ONLY");  
  dZ+=dl;
  gMC->Gspos("YMO2", 1, "YMOT", 0., 0., dZ, 0, "ONLY");  
  dZ+=dl;
//
//
// 2nd section: vacuum system 
//
  cpar0[0]=(kZvac7-kZvac4)/2;
  cpar0[1]=r2V;
  cpar0[2]=r2V+dVacuS;
  cpar0[3]=cpar0[1]+2.*cpar0[0]*TMath::Tan(kThetaOpenB);
  cpar0[4]=cpar0[2]+2.*cpar0[0]*TMath::Tan(kThetaOpenB);
  gMC->Gsvolu("YV21", "CONE", idtmed[kSteel+40], cpar0, 5);
//
// insulation
  cpar[0]=cpar0[0];
  cpar[1]=cpar0[1]+dTubeS;
  cpar[2]=cpar0[1]+dTubeS+kDInsuS;
  cpar[3]=cpar0[3]+dTubeS;
  cpar[4]=cpar0[3]+dTubeS+kDInsuS;
  gMC->Gsvolu("YI21", "CONE", idtmed[kInsulation+40], cpar, 5);
  gMC->Gspos("YI21", 1, "YV21", 0., 0., 0., 0, "ONLY"); 
  gMC->Gspos("YV21", 1, "YMO2", 0., 0., 0., 0, "ONLY"); 

//
// Third Section: Bellows and Flange 
//
  par3[0]  = 0.;
  par3[1]  = 360.;
  par3[2]  = 8.;
  dl=(kZvac9-kZvac7)/2.;
  
  par3[3]  = -dl;
  par3[4]  = r2+(kZvac7-kZvac3) * TMath::Tan(kThetaOpen2);
  par3[5]  = 30.;

  par3[6]  = -dl+kDr21;
  par3[7]  = par3[4]+kDr21;
  par3[8]  = 30.;

  par3[9]  = par3[6]+kDB2;
  par3[10] = par3[7];
  par3[11] = 30.;

  par3[12] = par3[9]+kDr22;
  par3[13] = par3[10]+kDr22;
  par3[14] = 30.;

  par3[15] = par3[12]+kDF2;
  par3[16] = par3[13];
  par3[17] = 30.;

  par3[18] = par3[15]+kDr22;
  par3[19] = par3[16]-kDr22;
  par3[20] = 30.;

  par3[21] = par3[18]+kDB2;
  par3[22] = par3[19];
  par3[23] = 30.;

  par3[24] = par3[21]+kDr23;
  par3[25] = par3[22];
  par3[26] = 30.;
//
  rBox=par3[22]-0.1;
  Float_t r3=par3[25];
  
  gMC->Gsvolu("YGO3", "PCON", idtmed[iHeavy+40], par3, 27);

  for (i=4; i<26; i+=3) par3[i]  = 0;

  gMC->Gsvolu("YMO3", "PCON", idtmed[kVacuum+40], par3, 27);
  gMC->Gspos("YGO3", 1, "YMO3", 0., 0., 0., 0, "ONLY");  

//
// Steel envelope
  tpar[0]=26;
  tpar[1]=30;
  tpar[2]=dl;
  gMC->Gsvolu("YS31", "TUBE", idtmed[kSteel], tpar, 3);
  gMC->Gspos("YS31", 1, "YGO3", 0., 0., 0., 0, "ONLY");  
  dZ+=dl;
  gMC->Gspos("YMO3", 1, "YMOT", 0., 0., dZ, 0, "ONLY");  
  dZ+=dl;

//
// 3rd section: vacuum system
//
//
// Bellow2
//
  tpar[0]=kRB2;
  tpar[1]=kRB2+kHB2;
  tpar[2]=kEB2/2.;
  gMC->Gsvolu("YB21", "TUBE", idtmed[kSteel+40], tpar, 3);
  dl1=tpar[2];
  
  tpar[0]=kRB2+kHB2-kEB2;
  tpar[1]=kRB2+kHB2;
  tpar[2]=(kLB2/2.-2.*kEB2)/2.;
  gMC->Gsvolu("YB22", "TUBE", idtmed[kSteel+40], tpar, 3);
  dl2=tpar[2];

  tpar[0]=kRB2-kEB2;
  tpar[1]=kRB2;
  tpar[2]=kLB2/8.;
  gMC->Gsvolu("YB23", "TUBE", idtmed[kSteel+40], tpar, 3);
  dl3=tpar[2];


  tpar[0]=0;
  tpar[1]=kRB2+kHB2;
  tpar[2]=kLB2/2.;
  gMC->Gsvolu("YBU2", "TUBE", idtmed[kVacuum+40], tpar, 3);

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
  tpar[1]=kRB2+kHB2;
  tpar[2]=7.*kLB2/2.;
  gMC->Gsvolu("YBM2", "TUBE", idtmed[kVacuum+40], tpar, 3);
  dz=-tpar[2]+kLB2/2.;

  for (i=0; i<7; i++) {
    gMC->Gspos("YBU2", i+1 , "YBM2", 0., 0.,dz , 0, "ONLY"); 
    dz+=kLB2;
  }

  dz=-dl+kDr21+tpar[2];
  gMC->Gspos("YBM2", 1, "YMO3", 0., 0., dz, 0, "ONLY"); 

  dz=dl-kDr23-tpar[2];
  gMC->Gspos("YBM2", 2, "YMO3", 0., 0., dz, 0, "ONLY"); 

//
// Flange

  tpar[0]=0;
  tpar[1]=kRF2;
  tpar[2]=kDF2/2.;
  gMC->Gsvolu("YFM2", "TUBE", idtmed[kVacuum+40], tpar, 3);

  tpar[0]=kRF2-2.;
  tpar[1]=kRF2;
  tpar[2]=kDF2/2.;
  gMC->Gsvolu("YF21", "TUBE", idtmed[kSteel+40], tpar, 3);
  gMC->Gspos("YF21", 1, "YFM2", 0., 0., 0., 0, "ONLY"); 

  tpar[0]=kRB2;
  tpar[1]=kRF2-2.;
  tpar[2]=kDFlange/2.;
  gMC->Gsvolu("YF22", "TUBE", idtmed[kSteel+40], tpar, 3);
  dz=-kDF2/2.+tpar[2];
  gMC->Gspos("YF22", 1, "YFM2", 0., 0., dz, 0, "ONLY"); 
  dz= kDF2/2.-tpar[2];
  gMC->Gspos("YF22", 2, "YFM2", 0., 0., dz, 0, "ONLY"); 

  dz=kDr21/2.-kDr23/2.;
  gMC->Gspos("YFM2", 2, "YMO3", 0., 0., dz, 0, "ONLY"); 


//
// pipe between flange and bellows
  tpar[0]=kRB2-dTubeS;
  tpar[1]=kRB2;
  tpar[2]=2.*(kDB2+kDr22-7.*kLB2)/4.;
  gMC->Gsvolu("YPF2", "TUBE", idtmed[kSteel+40], tpar, 3);
  dz=kDr21/2.-kDr23/2.-kDF2/2.-tpar[2];
  gMC->Gspos("YPF2", 1, "YMO3", 0., 0., dz, 0, "ONLY"); 
  dz=kDr21/2.-kDr23/2.+kDF2/2.+tpar[2];
  gMC->Gspos("YPF2", 2, "YMO3", 0., 0., dz, 0, "ONLY"); 

  Float_t dHorZ=20.;
  
//
// 4th section: rear shield and closing cone
//
  par4[0]  = 0.;
  par4[1]  = 360.;
  par4[2]  = 7.;
  dl=(kZvac12-kZvac9)/2.;
  
  par4[3]  = -dl;
  par4[4]  = r3;
  par4[5]  = 30.;

  par4[6]  = -dl+dHorZ;
  par4[7]  = r3;
  par4[8]  = 30.;

  par4[9]  = -dl+(kZvac10-kZvac9);
  par4[10]  = r3+(kZvac10-kZvac9-dHorZ) * TMath::Tan(kThetaOpen3);
  par4[11]  = 30.;

  par4[12]  = par4[9];
  par4[13] = par4[10];
  par4[14] = kR42;

  par4[15] = -dl+(kZvac11-kZvac9);
  par4[16] = r3+(kZvac11-kZvac9-dHorZ) * TMath::Tan(kThetaOpen3);
  par4[17] = kR42;

  par4[18] = par4[15];
  par4[19] = par4[16];
  par4[20] = kR43;

  par4[21] = -dl+(kZvac12-kZvac9);
  par4[22] = kRVacu+dVacuS;
  par4[23] = kR43;

  gMC->Gsvolu("YGO4", "PCON", idtmed[iHeavy+40], par4, 24);

//  parPb[0]  = (kZvac12-kZvac10)/2.;
//  parPb[1]  = parPb[3];
//  parPb[2]  = 31.;
//  parPb[3]  = parPb[1]+2.*parPb[0]*TMath::Tan(kThetaOpenPb);
//  parPb[4]  = 31.;
//  gMC->Gsvolu("YXO5", "CONE", idtmed[kPb], parPb, 5);
//  gMC->Gspos("YXO5", 1, "YGO4", 0., 0., -dl+(kZvac10-kZvac9)+parPb[0], 0, "ONLY");  

  for (i=4; i<23; i+=3) par4[i]  = 0;

  gMC->Gsvolu("YMO4", "PCON", idtmed[kVacuum+40], par4, 24);
  gMC->Gspos("YGO4", 1, "YMO4", 0., 0., 0., 0, "ONLY");  



  dZ+=dl;
  gMC->Gspos("YMO4", 1, "YMOT", 0., 0., dZ, 0, "ONLY");  
  dZ+=dl;
//
// Closing concrete cone 
//
  cpar[0]=(kZvac12-kZvac11)/2.;
  cpar[1] = r3+(kZvac11-kZvac9-dHorZ) * TMath::Tan(kThetaOpen3);
  cpar[2] = cpar[1]+0.001;
  cpar[3] = kRVacu+dVacuS;
  cpar[4] = cpar[2];
  gMC->Gsvolu("YCC4", "CONE", idtmed[kConcrete+40], cpar, 5);
  dz=dl-cpar[0];
  gMC->Gspos("YCC4", 1, "YGO4", 0., 0., dz, 0, "ONLY");  
//
// Steel envelope
//
  dz=-dl;
  tpar[0]=26.;
  tpar[1]=30.;
  tpar[2]=(kZvac10-kZvac9)/2.;
  gMC->Gsvolu("YS41", "TUBE", idtmed[kSteel], tpar, 3);
  dz+=tpar[2];
  gMC->Gspos("YS41", 1, "YGO4", 0., 0., dz, 0, "ONLY");  
  dz+=tpar[2];

  tpar[0]=kR41-kDRSteel2;
  tpar[1]=kR41;
  tpar[2]=(kZvac11-kZvac10)/2.;
  gMC->Gsvolu("YS43", "TUBE", idtmed[kPb], tpar, 3);
  dz+=tpar[2];
  gMC->Gspos("YS43", 1, "YGO4", 0., 0., dz, 0, "ONLY");  
//
// rear lead shield
//
  tpar[0]=kR41;
  tpar[1]=kR42;
  tpar[2]=(kZvac11-kZvac10)/2.;
  gMC->Gsvolu("YPBI", "TUBE", idtmed[kPb+40], tpar, 3);
  dz-=0;
  gMC->Gspos("YPBI", 1, "YGO4", 0., 0., dz, 0, "ONLY"); 

  tpar[0]=kR42-5;
  tpar[1]=kR42;
  tpar[2]=(kZvac11-kZvac10)/2.;
  gMC->Gsvolu("YPBO", "TUBE", idtmed[kPb], tpar, 3);
  gMC->Gspos("YPBO", 1, "YPBI", 0., 0., 0., 0, "ONLY"); 
  
//
// rear Fe shield
//

  tpar[0]=31.;
  tpar[1]=kR43;
  tpar[2]=(kZvac12-kZvac11)/2.;
  gMC->Gsvolu("YFEI", "TUBE", idtmed[kFe+40], tpar, 3);
  dz=dl-tpar[2];
  gMC->Gspos("YFEI", 1, "YGO4", 0., 0., dz, 0, "ONLY"); 

  tpar[0]=31.;
  tpar[1]=kR43;
  tpar[2]=2.5;
  gMC->Gsvolu("YFEO", "TUBE", idtmed[kFe], tpar, 3);
  dz=-(kZvac12-kZvac11)/2.+tpar[2];
  gMC->Gspos("YFEO", 1, "YFEI", 0., 0., dz, 0, "ONLY"); 
//
// Magnet element 
//
  tpar[0]= 0.;
  tpar[1]=40.;
  tpar[2]=85.;
  gMC->Gsvolu("YAEM", "TUBE", idtmed[kAir], tpar, 3);
  tpar[0]=17.6/2.;
  tpar[1]=40.;
  tpar[2]=85.;
  gMC->Gsvolu("YFEM", "TUBE", idtmed[kFe], tpar, 3);
  gMC->Gspos("YFEM", 1, "YAEM", 0., 0., 0., 0, "ONLY"); 

//
  dz=1921.6 + tpar[2];
  gMC->Gspos("YAEM", 1, "ALIC", 0., 0.,  dz, 0, "ONLY"); 

// 
//
// 4th section: vacuum system 
//
// up to closing cone
  
  Float_t r3V=r3-kDr23+dVacuS-1.6;

  cpar0[0]=(kZvac11-kZvac9)/2;
  cpar0[1]=r3V-dVacuS;
  cpar0[2]=r3V;
  cpar0[3]=cpar0[1]+2.*cpar0[0]*TMath::Tan(kThetaOpen3);
  cpar0[4]=cpar0[2]+2.*cpar0[0]*TMath::Tan(kThetaOpen3);
  gMC->Gsvolu("YV31", "CONE", idtmed[kSteel+40], cpar0, 5);
//
// insulation
  cpar[0]=cpar0[0];
  cpar[1]=cpar0[1]+dTubeS;
  cpar[2]=cpar0[1]+dTubeS+kDInsuS;
  cpar[3]=cpar0[3]+dTubeS;
  cpar[4]=cpar0[3]+dTubeS+kDInsuS;
  gMC->Gsvolu("YI31", "CONE", idtmed[kInsulation+40], cpar, 5);
  gMC->Gspos("YI31", 1, "YV31", 0., 0., 0., 0, "ONLY"); 
  dz=-dl+cpar[0];
  gMC->Gspos("YV31", 1, "YMO4", 0., 0., dz, 0, "ONLY"); 
//
// closing cone
  cpar0[0]=(kZvac12-kZvac11)/2;
  cpar0[1]=r3V-dVacuS+(kZvac11-kZvac9)*TMath::Tan(kThetaOpen3);
  cpar0[2]=r3V       +(kZvac11-kZvac9)*TMath::Tan(kThetaOpen3);
  cpar0[3]=kRVacu;
  cpar0[4]=kRVacu+dTubeS+kDInsuS+kDProtS+kDFreeS;
  gMC->Gsvolu("YV32", "CONE", idtmed[kSteel+40], cpar0, 5);
//
// insulation
  cpar[0]=cpar0[0];
  cpar[1]=cpar0[1]+dTubeS;
  cpar[2]=cpar0[1]+dTubeS+kDInsuS;
  cpar[3]=cpar0[3]+dTubeS;
  cpar[4]=cpar0[3]+dTubeS+kDInsuS;
  gMC->Gsvolu("YI32", "CONE", idtmed[kInsulation+40], cpar, 5);
  gMC->Gspos("YI32", 1, "YV32", 0., 0., 0., 0, "ONLY"); 
//
// clearance
//  cpar[1]=cpar0[2]-kDProtS-kDFreeS;
//  cpar[2]=cpar0[2]-kDProtS;
//  cpar[3]=cpar0[4]-kDProtS-kDFreeS;
//  cpar[4]=cpar0[4]-kDProtS;
//  gMC->Gsvolu("YP32", "CONE", idtmed[kVacuum+40], cpar, 5);
//  gMC->Gspos("YP32", 1, "YV32", 0., 0., 0., 0, "ONLY"); 
  
  dz=dl-cpar[0];
  gMC->Gspos("YV32", 1, "YMO4", 0., 0., dz, 0, "ONLY"); 
//
//
// MUON trigger wall
//  
  tpar[0] = 50.;
  tpar[1] = 310.;
  tpar[2] = (kZFilterOut - kZFilterIn) / 2.;
  gMC->Gsvolu("YFIM", "TUBE", idtmed[kFe+40], tpar, 3);
  dz = (kZFilterIn + kZFilterOut) / 2.;
  tpar[2] -= 10.;
  gMC->Gsvolu("YFII","TUBE", idtmed[kFe], tpar, 3);
  gMC->Gspos("YFII", 1, "YFIM", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("YFIM", 1, "ALIC", 0., 0., dz, 0, "ONLY");
//
// Shielding close to chamber
//
//
  cpar[0]=(kZch11-kZRear)/2.;
  cpar[1]=kR11;
  cpar[2]=kZRear*TMath::Tan(kAccMin);
  cpar[3]=kR11;
  cpar[4]=(kZRear+2.*cpar[0])*TMath::Tan(kAccMin);
  gMC->Gsvolu("YCS1", "CONE", idtmed[kNiCuW], cpar, 5);
  dz=-(kZvac12-zstart)/2.+(kZRear-zstart)+cpar[0];
  gMC->Gspos("YCS1", 1, "YMOT", 0., 0., dz, 0, "ONLY");

  cpar[0]=(kZvac4-kZch12)/2.;
  cpar[1]=kR11;
  cpar[2]=kZch12*TMath::Tan(kAccMin);
  cpar[3]=kR11;
  cpar[4]=(kZch12+2.*cpar[0])*TMath::Tan(kAccMin);
  gMC->Gsvolu("YCS3", "CONE", idtmed[kNiCuW], cpar, 5);
  dz=-(kZvac12-zstart)/2.+(kZch12-zstart)+cpar[0];
  gMC->Gspos("YCS3", 1, "YMOT", 0., 0., dz, 0, "ONLY");


// Recess station 1

  cpar[0]=(kZch12-kZch11)/2.;
  cpar[1]=kR11;
  cpar[2]=18.;
  cpar[3]=kR11;
  cpar[4]=17.9;
  gMC->Gsvolu("YCS2", "CONE", idtmed[kAir], cpar, 5);
  dz=-(kZvac12-zstart)/2.+(kZch11-zstart)+cpar[0];
  gMC->Gspos("YCS2", 1, "YMOT", 0., 0., dz, 0, "ONLY");

  Float_t ptubs[5];
  ptubs[0] = kR11;
  ptubs[1] = 17.9;
  ptubs[2] =   0.;
// phi_min, phi_max
  ptubs[3] =   0.;
  ptubs[4] =  90.;  
  gMC->Gsvolu("YCR0", "TUBS", idtmed[kNiCuW], ptubs, 0);
  Int_t idrotm[1799];
  
  AliMatrix(idrotm[1701],90.,   0., 90.,  90., 0., 0.);
  AliMatrix(idrotm[1702],90.,  90., 90., 180., 0., 0.);
  AliMatrix(idrotm[1703],90., 180., 90., 270., 0., 0.); 
  AliMatrix(idrotm[1704],90., 270., 90.,   0., 0., 0.); 
  //  Int_t ipos;
  
  dz=-cpar[0];
// 1.
  ptubs[2]=6.5/2.;
  dz+=ptubs[2];
  gMC->Gsposp("YCR0", 1, "YCS2", 0., 0., dz, idrotm[1701], "ONLY", ptubs, 5);
  gMC->Gsposp("YCR0", 2, "YCS2", 0., 0., dz, idrotm[1703], "ONLY", ptubs, 5);
  dz+=ptubs[2];
  dz+=1.5;
// 2.
  ptubs[2]=5.0/2.;
  dz+=ptubs[2];
  gMC->Gsposp("YCR0", 3, "YCS2", 0., 0., dz, idrotm[1702], "ONLY", ptubs, 5);
  gMC->Gsposp("YCR0", 4, "YCS2", 0., 0., dz, idrotm[1704], "ONLY", ptubs, 5);
  dz+=ptubs[2];
  dz+=1.5;
// 3. 
  ptubs[2]=5.0/2.;
  dz+=ptubs[2];
  gMC->Gsposp("YCR0", 5, "YCS2", 0., 0., dz, idrotm[1701], "ONLY", ptubs, 5);
  gMC->Gsposp("YCR0", 6, "YCS2", 0., 0., dz, idrotm[1703], "ONLY", ptubs, 5);
  dz+=ptubs[2];
  dz+=1.5;
// 4. 
  ptubs[2]=6.5/2.;
  dz+=ptubs[2];
  gMC->Gsposp("YCR0", 7, "YCS2", 0., 0., dz, idrotm[1702], "ONLY", ptubs, 5);
  gMC->Gsposp("YCR0", 8, "YCS2", 0., 0., dz, idrotm[1704], "ONLY", ptubs, 5);
  dz+=ptubs[2];
  dz+=1.5;


  
  cpar[0]=(kZch21-kZvac4)/2.;
  cpar[1]=kR21;
  cpar[2]=kZvac4*TMath::Tan(kAccMin);
  cpar[3]=kR21;
  cpar[4]=(kZvac4+2.*cpar[0])*TMath::Tan(kAccMin);
  gMC->Gsvolu("YCS4", "CONE", idtmed[kNiCuW], cpar, 5);
  dz=-(kZvac12-zstart)/2.+(kZvac4-zstart)+cpar[0];
  gMC->Gspos("YCS4", 1, "YMOT", 0., 0., dz, 0, "ONLY");

  cpar[0]=(kZvac6-kZch22)/2.;
  cpar[1]=kR21;
  cpar[2]=kZch22*TMath::Tan(kAccMin);
  cpar[3]=kR21;
  cpar[4]=(kZch22+2.*cpar[0])*TMath::Tan(kAccMin);
  gMC->Gsvolu("YCS6", "CONE", idtmed[kNiCuW], cpar, 5);
  dz=-(kZvac12-zstart)/2.+(kZch22-zstart)+cpar[0];
  gMC->Gspos("YCS6", 1, "YMOT", 0., 0., dz, 0, "ONLY");
  
// Recess station 2
 
  cpar[0]=(kZch22-kZch21)/2.;
  cpar[1]=kR21;
  cpar[2]=23.;
  cpar[3]=kR21;
  cpar[4]=23.;
  gMC->Gsvolu("YCS5", "CONE", idtmed[kAir], cpar, 5);
  dz=-(kZvac12-zstart)/2.+(kZch21-zstart)+cpar[0];
  gMC->Gspos("YCS5", 1, "YMOT", 0., 0., dz, 0, "ONLY");

  ptubs[0] = kR21;
  ptubs[1] = 23;
  ptubs[2] =   0.;
  ptubs[3] =   0.;
  ptubs[4] =  90.;  
  gMC->Gsvolu("YCR1", "TUBS", idtmed[kNiCuW], ptubs, 0);

  dz=-cpar[0];
// 1.
  ptubs[2]=7.5/2.;
  dz+=ptubs[2];
  gMC->Gsposp("YCR1", 1, "YCS5", 0., 0., dz, idrotm[1701], "ONLY", ptubs, 5);
  gMC->Gsposp("YCR1", 2, "YCS5", 0., 0., dz, idrotm[1703], "ONLY", ptubs, 5);
  dz+=ptubs[2];
  dz+=1.5;
// 2.
  ptubs[2]=6.0/2.;
  dz+=ptubs[2];
  gMC->Gsposp("YCR1", 3, "YCS5", 0., 0., dz, idrotm[1702], "ONLY", ptubs, 5);
  gMC->Gsposp("YCR1", 4, "YCS5", 0., 0., dz, idrotm[1704], "ONLY", ptubs, 5);
  dz+=ptubs[2];
  dz+=1.5;
// 3. 
  ptubs[2]=6.0/2.;
  dz+=ptubs[2];
  gMC->Gsposp("YCR1", 5, "YCS5", 0., 0., dz, idrotm[1701], "ONLY", ptubs, 5);
  gMC->Gsposp("YCR1", 6, "YCS5", 0., 0., dz, idrotm[1703], "ONLY", ptubs, 5);
  dz+=ptubs[2];
  dz+=1.5;
// 4. 
  ptubs[2]=7.5/2.;
  dz+=ptubs[2];
  gMC->Gsposp("YCR1", 7, "YCS5", 0., 0., dz, idrotm[1702], "ONLY", ptubs, 5);
  gMC->Gsposp("YCR1", 8, "YCS5", 0., 0., dz, idrotm[1704], "ONLY", ptubs, 5);
  dz+=ptubs[2];
  dz+=1.5;

//
// Outer Pb Cone

  if (fPbCone) {
      dl = (kZvac10-kZch32)/2.;
      dz = dl+kZch32;
      
      par0[0]  = 0.;
      par0[1]  = 360.;
      par0[2]  = 10.;

      par0[ 3]  = -dl;
      par0[ 4]  = 30.;
      par0[ 5]  = 30.+(kZch32-kZConeE)*TMath::Tan(kThetaOpenPbO);

//    4th station
      par0[ 6]  = -dz + kZch41;
      par0[ 7]  = 30.;
      par0[ 8]  = 30.+(kZch41-kZConeE)*TMath::Tan(kThetaOpenPbO);

      par0[ 9]  = -dz + kZch41;
      par0[10]  = 30.;
      par0[11]  = 37.5;  
                                          // recess erice2000
      par0[12]  = -dz + kZch42;
      par0[13]  = 30.;
      par0[14]  = par0[11];

      par0[15]  = -dz + kZch42;
      par0[16]  = 30.;
      par0[17]  = 30.+(kZch42-kZConeE)*TMath::Tan(kThetaOpenPbO);

//    5th station
      par0[18]  = -dz + kZch51;
      par0[19]  = 30.;
      par0[20]  = 30.+(kZch51-kZConeE)*TMath::Tan(kThetaOpenPbO);

      par0[21]  = -dz + kZch51;
      par0[22]  = 30.;
      par0[23]  = 37.5;  // recess erice2000

      par0[24]  = -dz + kZch52;
      par0[25]  = 30.;
      par0[26]  = par0[23];

      par0[27]  = -dz + kZch52;
      par0[28]  = 30.;
      par0[29]  = 30.+(kZch52-kZConeE)*TMath::Tan(kThetaOpenPbO);
// end of cone
      par0[30]  = +dl;
      par0[31]  = 30.;
      par0[32]  = par0[29];
//
      gMC->Gsvolu("YOPB", "PCON", idtmed[kPb], par0, 33);
      dz = -(kZvac12-zstart)/2. + (kZch32-zstart) + dl;
      gMC->Gspos("YOPB", 1, "YMOT", 0., 0., dz, 0, "ONLY");
  }
}

void AliSHILv0::Init()
{
  //
  // Initialise the muon shield after it has been built
  //
  Int_t i;
  //
  
  if(AliLog::GetGlobalDebugLevel()>0) {
      printf("\n%s: ",ClassName());
      for(i=0;i<35;i++) printf("*");
      printf(" SHILv0_INIT ");
      for(i=0;i<35;i++) printf("*");
      printf("\n%s: ",ClassName());
      //
      // Here the SHIL initialisation code (if any!)
      for(i=0;i<80;i++) printf("*");
      printf("\n");
  }
}
