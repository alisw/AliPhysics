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

//-------------------------------------------------------------------------
// MUON shielding class
// Default version
// Author: A.Morsch
//-------------------------------------------------------------------------

#include <TVirtualMC.h>
#include <TArrayI.h>
#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>
#include <TGeoBBox.h>
#include <TGeoPgon.h>
#include <TGeoTorus.h>

#include "AliSHILv2.h"
#include "AliConst.h"
#include "AliLog.h"

ClassImp(AliSHILv2)

 
//_____________________________________________________________________________

AliSHILv2::AliSHILv2():
    fPbCone(kTRUE),
    fWriteGeometry(kFALSE)
{
  //
  // Default constructor for muon shield
  //
}
 
//_____________________________________________________________________________
AliSHILv2::AliSHILv2(const char *name, const char *title)
    : AliSHIL(name,title),
      fPbCone(kTRUE),
      fWriteGeometry(kFALSE)
{
  //
  // Standard constructor for muon shield
  //
  // Pb  cone not yet compatible with muon chamber inner radii
  // Switched off by default
}
 
//_____________________________________________________________________________
void AliSHILv2::CreateGeometry()
{
  // 
  // Build muon shield geometry
  //
  //
  //Begin_Html
  /*
    <img src="picts/AliSHILv2.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliSHILv2Tree.gif">
  */
  //End_Html

    Float_t cpar[5], cpar0[5], tpar[3], par1[100], pars1[100], par2[100], par3[100], 
	par4[24], par0[100];
    Float_t dz, dZ;
    
    Int_t *idtmed = fIdtmed->GetArray()-1699;

    Int_t idrotm[1799];

#include "ABSOSHILConst.h"
#include "SHILConst2.h"
    
    enum {kC=1705, kAl=1708, kFe=1709, kCu=1710, kW=1711, kPb=1712,
	  kNiCuW=1720, kVacuum=1715, kAir=1714, kConcrete=1716,
	  kPolyCH2=1717, kSteel=1718, kInsulation=1713, kAirMuon = 1774};	
    Int_t i;
    
//
// Material of the rear part of the shield
  Int_t iHeavy = kNiCuW;
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
  const Float_t kzLength = dl;
  
  dz=zstart+dl;
//
// start
  par0[3]  = -dl;
  par0[4]  = 0.;
  par0[5]  = zstart * TMath::Tan(kAccMin);
// recess station 1
  par0[6]  = -dl  - zstart + kZch11;
  par0[7]  = 0.;
  par0[8]  = 18.2;

  par0[9]   = par0[6];
  par0[10]  = 0.;
  par0[11]  = kR11;

  par0[12]  = -dl - zstart + kZch12;
  par0[13]  = 0.;
  par0[14]  = kR11;

  par0[15]  = par0[12];
  par0[16]  = 0.;
  par0[17]  = 19.5;
// recess station 2
  par0[18]  = -dz+kZch21;
  par0[19]  = 0.;
  par0[20]  = kZch21 * TMath::Tan(kAccMin);

  par0[21]  = -dz+kZch21;
  par0[22] = 0.;
  par0[23] = kR21;

  par0[24]  = -dz+kZch22;
  par0[25] = 0.;
  par0[26] = kR21;

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
  par0[41] = 28.8;

  par0[42] = -dz+kZch32;
  par0[43] = 0.;
  par0[44] = 28.8;
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
  par0[53] = 36.9;

  par0[54] = -dz+kZch42;
  par0[55] = 0.;
  par0[56] = 36.9;

  par0[57] = -dz+kZch42;
  par0[58] = 0.;
  par0[59] = 30.+(kZch42-kZConeE)*TMath::Tan(kThetaOpenPbO);

// recess station 5

  par0[60] = -dz+kZch51;
  par0[61] = 0.;
  par0[62] = 30.+(kZch51-kZConeE)*TMath::Tan(kThetaOpenPbO);

  par0[63] = -dz+kZch51;
  par0[64] = 0.;
  par0[65] = 36.9;

  par0[66] = -dz+kZch52;
  par0[67] = 0.;
  par0[68] = 36.9;

  par0[69] = -dz+kZch52;
  par0[70] = 0.;
  par0[71] =  30.+(kZch52+4.-kZConeE)*TMath::Tan(kThetaOpenPbO);

// end of cone

  par0[72] = -dz+kZvac10;
  par0[73] = 0.;
  par0[74] = par0[71];

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
  AliMatrix(idrotm[1705], 270., 0., 90., 90., 180., 0.);

  dZ=-dl;

//
// Dimuon arm mother volumes (YOUT1, YOUT2)
//

  // Dipole parameters
  Float_t z01 =  -724.45;
  // Float_t z02 =  -814.30;
  Float_t z05 = -1235.55;

// Before dipole
//

  Float_t zpos  = -zstart - kzLength;
  Float_t shift_after_absorber = 35;
  Float_t delta = 1e-06;
  Float_t rst1  = 120;
  Float_t rst2  = 150;
  Float_t rst3  = 252;
  Float_t rst4  = 252;
  Float_t rst5  = 304;
  //Float_t rst6  = 430.;
  Float_t rst7  = 460.;

  // Float_t zstart2 = zpos + zstart;

  par0[0]  = 0.;
  par0[1]  = 360.;
  par0[2]  = 11.;

//
// start
  // par0[3]  = zpos - ( -dl );
       // start only after absorber
  // z = -503.00     
  par0[3]  = zpos - ( -dl + shift_after_absorber ) + delta;
  par0[4]  = ( zstart + shift_after_absorber ) * TMath::Tan(kAccMin) + delta;
  par0[5]  = rst1;

// recess station 1
  // z = -517.70     
  par0[6]  = zpos - ( -dl  - zstart + kZch11 + delta );
  par0[7]  = 18.2 + delta;
  par0[8]  = rst1;

  // z = -517.70     
  par0[9]   = par0[6];
  par0[10]  = kR11 + delta;
  par0[11]  = rst1;

  // z = -553.70     
  par0[12]  = zpos - ( -dl - zstart + kZch12 - delta );
  par0[13]  = kR11 + delta;
  par0[14]  = rst1;

  // z = -553.70     
  par0[15]  = par0[12];
  par0[16]  = 19.5 + delta;
  par0[17]  = rst1;

// recess station 2
  // z = -661.30     
  par0[18]  = zpos - ( -dz+kZch21 + delta );
  par0[19]  = kZch21 * TMath::Tan(kAccMin) + delta;
  par0[20]  = rst2;

  // z = -661.30     
  par0[21]  = par0[18];
  par0[22] = kR21 + delta;
  par0[23] = rst2;

  // z = -709.90     
  par0[24]  = zpos - ( -dz+kZch22 - delta );
  par0[25] = kR21 + delta;
  par0[26] = rst2;

  // z = -709.90     
  par0[27]  = par0[24];
  par0[28]  = kZch22 * TMath::Tan(kAccMin) + delta;
  par0[29]  = rst2;
//
  // z = -711.00     
  par0[30] = zpos - ( -dz+kZvac6 );
  par0[31] = kZvac6 * TMath::Tan(kAccMin) + delta;
  par0[32] = rst2;

  Float_t nextZ = zpos - ( -dz+kZConeE );
  Float_t nextRin = 30. + delta;
  Float_t nextRout = rst3;
  Float_t tgin  = ( nextRin - par0[31]) / (nextZ - par0[30]);
  Float_t tgout = ( nextRout - par0[32])/ (nextZ - par0[30]);

  // z = -724.45     
  par0[33] = z01;
  par0[34] = par0[31] + (z01 - par0[30]) * tgin;
  par0[35] = par0[32] + (z01 - par0[30]) * tgout;

  gMC->Gsvolu("YOUT1", "PCON", idtmed[kAirMuon], par0, 36);
  gMC->Gspos("YOUT1", 1, "ALIC", 0., 0., 0., 0, "ONLY");  


// 
// After dipole
// 

  par0[0]  = 0.;
  par0[1]  = 360.;
  par0[2]  = 14.;

  // z = -1235.55
  par0[3] = z05;
  par0[4] = nextRin - (z05 - nextZ) * TMath::Tan(kThetaOpenPbO) + delta ;
  par0[5] = rst4;
    
// recess station 4
  // z = -1259.90
  par0[6] = zpos - ( -dz+kZch41 + delta );
  par0[7] = 30.+(kZch41-kZConeE)*TMath::Tan(kThetaOpenPbO) + delta;
  par0[8] = rst4;

  // z = -1259.90
  par0[9] = par0[6];
  par0[10] = 36.9 + delta;
  par0[11] = rst4;

  // z = -1324.10
  par0[12] = zpos - ( -dz+kZch42 - delta );
  par0[13] = 36.9 + delta;
  par0[14] = rst4;

  // z = -1324.10
  par0[15] = par0[12];
  par0[16] = 30.+(kZch42-kZConeE)*TMath::Tan(kThetaOpenPbO) + delta;
  par0[17] = rst5;

// recess station 5

  // z = -1390.00
  par0[18] = zpos - ( -dz+kZch51 + delta );
  par0[19] = 30.+(kZch51-kZConeE)*TMath::Tan(kThetaOpenPbO) + delta;
  par0[20] = rst5;

  // z = -1390.00
  par0[21] = par0[18];
  par0[22] = 36.9 + delta;
  par0[23] = rst5;

  // z = -1454.20
  par0[24] = zpos - ( -dz+kZch52 - delta );
  par0[25] = 36.9 + delta;
  par0[26] = rst5;

  // z = -1454.20
  par0[27] = par0[24];
  par0[28] =  30.+(kZch52+4.-kZConeE)*TMath::Tan(kThetaOpenPbO) + delta;
  par0[29] = rst5;

// end of cone

  // z = -1466.00
  par0[30] = zpos - ( -dz+kZvac10 - delta );
  par0[31] = par0[28];
  par0[32] = rst7;

  // z = -1466.00
  par0[33] = par0[30];
  par0[34] = kR42 + delta;
  par0[35] = rst7;

  // z = -1800.00
  par0[36] = zpos - ( -dz+kZvac11 - delta );
  par0[37] = kR42 + delta;
  par0[38] = rst7;

  // z = -1800.00
  par0[39] = par0[36];
  par0[40] = kR43 + delta;
  par0[41] = rst7;

  // z = -1900.00
  par0[42] = zpos - ( -dz+kZvac12 );
  par0[43] = kR43 + delta;
  par0[44] = rst7;

  gMC->Gsvolu("YOUT2", "PCON", idtmed[kAirMuon], par0, 45);
  gMC->Gspos("YOUT2", 1, "ALIC", 0., 0., 0., 0, "ONLY");  

//
// First section: bellows below and behind front absorber 
// 
//
  par1[ 0]  = 0.;
  par1[ 1]  = 360.;
  par1[ 2]  = 15.;
  dl=(kZvac4-zstart)/2.;
  
  par1[ 3]  = -dl;
  par1[ 4]  = kRAbs+(zstart-kZOpen) * TMath::Tan(kThetaOpen1);
  par1[ 5]  = zstart * TMath::Tan(kAccMin);

  par1[ 6]  = -dl-zstart+kZch11;
  par1[ 7]  = par1[4] + (dRear1 + 19.)  * TMath::Tan(kThetaOpen1);
  par1[ 8]  = 18.2;

  par1[ 9]  = par1[6];
  par1[10]  = par1[7];
  par1[11]  = kR11;

  par1[12]  = -dl-zstart+kZch12;
  par1[13]  = par1[10] + 36. * TMath::Tan(kThetaOpen1);
  par1[14]  = kR11;

  par1[15]  = -dl+dRear1 + 50.7;
  par1[16]  = par1[13];
  par1[17]  = 19.5;

  par1[18]  = -dl+kZvac1-zstart;
  par1[19]  = par1[16] + (par1[18] - par1[15]) * TMath::Tan(kThetaOpen1);
  par1[20]  = (par1[18] +dl +zstart) * TMath::Tan(kAccMin);

  par1[21]  = -dl+kZvac1-zstart;
  par1[22]  = kRAbs+ (kZvac1-kZOpen) * TMath::Tan(kThetaOpen1);
  par1[23]  = (par1[21] +dl +zstart) * TMath::Tan(kAccMin);
  

  par1[24]  = par1[21]+kDr11/10.;
  par1[25]  = par1[22]+kDr11;
  par1[26]  = (par1[24] +dl +zstart) * TMath::Tan(kAccMin);

  par1[27]  = -dl+(kZvac1+kDr11/10.+kDB1-zstart);
  par1[28]  = par1[25];
  par1[29]  = (par1[27] +dl +zstart) * TMath::Tan(kAccMin);

  par1[30]  = par1[27]+kDr12;
  par1[31]  = par1[28]+kDr12;
  par1[32]  = (par1[30] +dl +zstart) * TMath::Tan(kAccMin);

  par1[33]  = par1[30]+kDF1;
  par1[34]  = par1[31];
  par1[35]  = (par1[33] +dl +zstart) * TMath::Tan(kAccMin);

  par1[36]  = par1[33]+kDr12;
  par1[37]  = par1[34]-kDr12; 
  par1[38]  = (par1[36] +dl +zstart) * TMath::Tan(kAccMin);

  par1[39] = par1[36]+kDB1;
  par1[40] = par1[37];
  par1[41] = (par1[39] +dl +zstart) * TMath::Tan(kAccMin);

  par1[42] = par1[39]+kDr13;
  par1[43] = par1[40]-kDr13;
  par1[44] = (par1[42] +dl +zstart) * TMath::Tan(kAccMin);

  par1[45] =  -dl+kZvac4-zstart;
  par1[46] = par1[43];
  par1[47] = (par1[45] +dl +zstart) * TMath::Tan(kAccMin);

  Float_t r2  = par1[46];
  Float_t rBox= par1[46]-0.1;

  gMC->Gsvolu("YGO1", "PCON", idtmed[kNiCuW+40], par1, 48);

  for (i=0; i<48; i++)  pars1[i]  = par1[i];
  for (i=4; i<47; i+=3) pars1[i]  = 0.;

  gMC->Gsvolu("YMO1", "PCON", idtmed[kVacuum+40], pars1, 48);
  gMC->Gspos("YGO1", 1, "YMO1", 0., 0., 0., 0, "ONLY");  
  dZ+=dl;
  gMC->Gspos("YMO1", 1, "YMOT", 0., 0., dZ, 0, "ONLY");  
  dZ+=dl;


  tpar[0]=kR21-0.6;
  tpar[1]=kR21;
  tpar[2]=(kZvac4-kZvac41)/2.;
  gMC->Gsvolu("YSE1", "TUBE", idtmed[kSteel], tpar, 3);
  dz=dl-tpar[2];
  gMC->Gspos("YSE1", 1, "YGO1", 0., 0., dz, 0, "ONLY");


  tpar[0]=kR11-0.6;
  tpar[1]=kR11;
  tpar[2]=(kZvac41-zstart-dRear1)/2.;
  gMC->Gsvolu("YSE2", "TUBE", idtmed[kSteel], tpar, 3);
  dz=dl-tpar[2]-(kZvac4-kZvac41);
  gMC->Gspos("YSE2", 1, "YGO1", 0., 0., dz, 0, "ONLY");

//
// 1st section: vacuum system
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
  tpar[2]=kLB1/2.;
  gMC->Gsvolu("YBU1", "TUBE", idtmed[kVacuum+40], tpar, 3);

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
  tpar[1]=kRB1+kHB1+0.5;
  tpar[2]=12.*kLB1/2.;
  gMC->Gsvolu("YBM1", "TUBE", idtmed[kVacuum+40], tpar, 3);
  Float_t bsize = tpar[2];
  tpar[0]=kRB1+kHB1;
  gMC->Gsvolu("YBI1", "TUBE", idtmed[kInsulation+40], tpar, 3);
  gMC->Gspos("YBI1", 2, "YBM1", 0., 0., 0., 0, "ONLY"); 

  dz=-bsize+kLB1/2.;

  for (i=0; i<12; i++) {
    gMC->Gspos("YBU1", i+1 , "YBM1", 0., 0., dz, 0, "ONLY"); 
    dz+=kLB1;
  }

  dz=-dl+(kZvac1-zstart)+kDr11/10.+bsize;
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
  dz=-dl+(kZvac3-zstart)-2.*kDr13-tpar[2];
  gMC->Gspos("YFM1", 2, "YMO1", 0., 0., dz, 0, "ONLY"); 

//
// pipe between flange and bellows
//
// Steel 
  tpar[0] = kRB1-dTubeS;
  tpar[1] = kRB1+0.6;
  tpar[2] = (kZvac3-kZvac1-2.*kDr13-kDr11/10.-kDF1-2.*bsize)/2.;
  gMC->Gsvolu("YPF1", "TUBE", idtmed[kSteel+40], tpar, 3);
// Insulation
  tpar[0]=kRB1;
  tpar[1]=kRB1+0.5;
  gMC->Gsvolu("YPS1", "TUBE", idtmed[kInsulation+40], tpar, 3);
  gMC->Gspos("YPS1", 1, "YPF1", 0., 0., 0., 0, "ONLY"); 
  dz=-dl+(kZvac1-zstart)+kDr11/10.+2.*bsize+tpar[2];
  gMC->Gspos("YPF1", 1, "YMO1", 0., 0., dz, 0, "ONLY"); 


// Pipe+Heating     1.5 mm 
// Heating Jacket   5.0 mm
// Protection       1.0 mm
// ========================
//                  7.5 mm
// pipe and heating jackets outside bellows
//
// left side
  cpar0[0]=(kZvac1+kDr11/10.-zstart)/2;
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
  
  cpar0[0] = (kZvac4-kZvac3+2.*kDr13)/2;
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
  par2[26] = 28.8;

  par2[27] = -dl+(kZch32-kZvac4);
  par2[28] = r2+(kZch32-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  par2[29] = 28.8;

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
  Float_t parPb[18];
  parPb[ 0]  = 0.;
  parPb[ 1]  = 360.;
  parPb[ 2]  = 5.;
  Float_t dlPb=(kZvac7-kZPb)/2.;
  
  parPb[ 3]  = -dlPb;
  parPb[ 4]  =  r2+(kZPb-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  parPb[ 5]  =  kZPb*TMath::Tan(kAccMin)-kDRSteel2;
  
  parPb[ 6]  = -dlPb+(kZConeE-kZPb);
  parPb[ 7]  =  r2+(kZConeE-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  parPb[ 8]  = 26.;

  parPb[ 9]  = -dlPb+(kZch32+4.-kZPb);
  parPb[10]  =  r2+(kZch32+4.-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  parPb[11]  = 26.;

  parPb[12]  = -dlPb+(kZch32+4.-kZPb);
  parPb[13]  =  r2+(kZch32+4.-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  parPb[14]  = 30.;
  
  parPb[15]  = dlPb;
  parPb[16]  =  r2+(kZvac7-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  parPb[17]  = 30.;

  gMC->Gsvolu("YXO2", "PCON", idtmed[kPb+40], parPb, 18);	  
  gMC->Gspos("YXO2", 1, "YGO2", 0., 0., (kZPb-kZvac4)/2., 0, "ONLY");  
//
// Concrete replacing Pb
//
  Float_t parCC[9];
  Float_t zCC1 = 1066.;
  Float_t zCC2 = 1188.;
  
  parCC[ 0]  = 0.;
  parCC[ 1]  = 360.;
  parCC[ 2]  = 2.;
  Float_t dlCC=(zCC2-zCC1)/2.;
  parCC[ 3]  = -dlCC;
  parCC[ 4]  =  r2+(zCC1-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  parCC[ 5]  =  30.;
  
  parCC[ 6]  =  dlCC;
  parCC[ 7]  =  r2+(zCC2-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  parCC[ 8]  = 30.;
  gMC->Gsvolu("YCO2", "PCON", idtmed[kSteel], parCC, 9);	  
//  gMC->Gspos("YCO2", 1, "YXO2", 0., 0., dlPb-dlCC-(kZvac7-zCC2), 0, "ONLY");  

  zCC1 = 751.75;
  zCC2 = kZConeE;
  dlCC=(zCC2-zCC1)/2.;
  parCC[ 3]  = -dlCC;
  parCC[ 4]  =  r2+(zCC1-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  parCC[ 5]  =  zCC1*TMath::Tan(kAccMin)-kDRSteel2;
  
  parCC[ 6]  =  dlCC;
  parCC[ 7]  =  r2+(zCC2-kZvac4-10.) * TMath::Tan(kThetaOpen2);
  parCC[ 8]  = 26.;
  
  gMC->Gsvolu("YCO1", "PCON", idtmed[kSteel], parCC, 9);	  
//  gMC->Gspos("YCO1", 1, "YXO2", 0., 0., dlPb-dlCC-(kZvac7-zCC2), 0, "ONLY");  
  
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
 
  parW[12]   = dlW;
  parW[13]  =  r2+(kZPb-kZvac4) * TMath::Tan(kThetaOpen2);
  parW[14]  = kZPb*TMath::Tan(kAccMin)-kDRSteel2;

  gMC->Gsvolu("YYO2", "PCON", idtmed[kNiCuW+40], parW, 15);	  
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
//  tpar[0]=26;
//  tpar[1]=30;
//  tpar[2]=dl;
//  gMC->Gsvolu("YS31", "TUBE", idtmed[kSteel], tpar, 3);
//  gMC->Gspos("YS31", 1, "YGO3", 0., 0., 0., 0, "ONLY");  
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


  for (i=4; i<23; i+=3) par4[i]  = 0;

  gMC->Gsvolu("YMO4", "PCON", idtmed[kVacuum+40], par4, 24);
  gMC->Gspos("YGO4", 1, "YMO4", 0., 0., 0., 0, "ONLY");  



  dZ+=dl;
  gMC->Gspos("YMO4", 1, "YMOT", 0., 0., dZ, 0, "ONLY");  
  dZ+=dl;
//
// Concrete replacing Pb
//
  zCC1 = 1316.;
  zCC2 = 1349.;
  
  parCC[ 0]  = 0.;
  parCC[ 1]  = 360.;
  parCC[ 2]  = 2.;
  dlCC=(zCC2-zCC1)/2.;
  parCC[ 3]  = -dlCC;
  parCC[ 4]  = r3+(zCC1-kZvac9-dHorZ) * TMath::Tan(kThetaOpen3);
  parCC[ 5]  =  30.;
  
  parCC[ 6]  =  dlCC;
  parCC[ 7]  =  r3+(zCC2-kZvac9-dHorZ) * TMath::Tan(kThetaOpen3);
  parCC[ 8]  = 30.;

  gMC->Gsvolu("YCO4", "PCON", idtmed[kSteel], parCC, 9);	  
//  gMC->Gspos("YCO4", 1, "YGO4", 0., 0., dl-dlCC-(kZvac12-zCC2), 0, "ONLY");  

  zCC1 = 1471.;
  zCC2 = 1591.;

  dlCC=(zCC2-zCC1)/2.;
  parCC[ 3]  = -dlCC;
  parCC[ 4]  = r3+(zCC1-kZvac9-dHorZ) * TMath::Tan(kThetaOpen3);
  parCC[ 5]  = kR41-kDRSteel2;
  
  parCC[ 6]  =  dlCC;
  parCC[ 7]  =  r3+(zCC2-kZvac9-dHorZ) * TMath::Tan(kThetaOpen3);
  parCC[ 8]  =  kR41-kDRSteel2;

  gMC->Gsvolu("YCO5", "PCON", idtmed[kSteel], parCC, 9);	  
//  gMC->Gspos("YCO5", 1, "YGO4", 0., 0., dl-dlCC-(kZvac12-zCC2), 0, "ONLY");  

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
//  gMC->Gspos("YS41", 1, "YGO4", 0., 0., dz, 0, "ONLY");  
  dz+=tpar[2];

  tpar[0]=kR41-kDRSteel2;
  tpar[1]=kR41;
  tpar[2]=(kZvac11-kZvac10)/2.;
  gMC->Gsvolu("YS43", "TUBE", idtmed[kPb+40], tpar, 3);
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

  tpar[2]=(zCC2-zCC1)/2.;
  gMC->Gsvolu("YCO6", "TUBE", idtmed[kSteel+40], tpar, 3);
//  gMC->Gspos("YCO6", 1, "YPBI", 0., 0., -(kZvac11-kZvac10)/2.+tpar[2], 0, "ONLY"); 


  tpar[0]=kR42-5;
  tpar[1]=kR42;
  tpar[2]=(kZvac11-kZvac10)/2.;
  gMC->Gsvolu("YPBO", "TUBE", idtmed[kPb+40], tpar, 3);
  gMC->Gspos("YPBO", 1, "YPBI", 0., 0., 0., 0, "ONLY"); 

  tpar[2]=(zCC2-zCC1)/2.;
  gMC->Gsvolu("YCO7", "TUBE", idtmed[kSteel], tpar, 3);
//  gMC->Gspos("YCO7", 1, "YPBO", 0., 0., -(kZvac11-kZvac10)/2.+tpar[2], 0, "ONLY"); 
  
//
// rear Fe shield
//

  tpar[0]=31.;
  tpar[1]=kR43;
  tpar[2]=(kZvac12-kZvac11)/2.;
  gMC->Gsvolu("YFEI", "TUBE", idtmed[kSteel+40], tpar, 3);
  dz=dl-tpar[2];
  gMC->Gspos("YFEI", 1, "YGO4", 0., 0., dz, 0, "ONLY"); 

  tpar[0]=31.;
  tpar[1]=kR43;
  tpar[2]=2.5;
  gMC->Gsvolu("YFEO", "TUBE", idtmed[kSteel], tpar, 3);
  dz=-(kZvac12-kZvac11)/2.+tpar[2];
  gMC->Gspos("YFEO", 1, "YFEI", 0., 0., dz, 0, "ONLY"); 


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
  
  dz=dl-cpar[0];
  gMC->Gspos("YV32", 1, "YMO4", 0., 0., dz, 0, "ONLY"); 

///////////////////////////////////
//    Muon Filter                //
//    Drawing ALIP2A__0105       //
///////////////////////////////////
  TGeoBBox*   shMuonFilterO1  = new TGeoBBox(550./2., 620./2., 120./2.);
  shMuonFilterO1->SetName("FilterO1");
  TGeoTube*   shMuonFilterI1  = new TGeoTube(0., 50., 121./2.);
  shMuonFilterI1->SetName("FilterI1");
  TGeoCompositeShape* shMuonFilterM = new TGeoCompositeShape("YMuonFilterM", "FilterO1-FilterI1");
  TGeoVolume* voMuonFilterM = new TGeoVolume("YMuonFilterM", shMuonFilterM,  gGeoManager->GetMedium("SHIL_ST_C0"));
  
  TGeoBBox*   shMuonFilterO2  = new TGeoBBox(550./2., 620./2., 110./2.);
  shMuonFilterO2->SetName("FilterO2");
  TGeoTube*   shMuonFilterI2  = new TGeoTube(0., 50., 111./2.);
  shMuonFilterI2->SetName("FilterI2");
  TGeoCompositeShape* shMuonFilterI = new TGeoCompositeShape("YMuonFilterI", "FilterO2-FilterI2");
  TGeoVolume* voMuonFilterI = new TGeoVolume("YMuonFilterI", shMuonFilterI,  gGeoManager->GetMedium("SHIL_ST_C3"));
  voMuonFilterI->SetName("YMuonFilterI");
  voMuonFilterM->AddNode(voMuonFilterI, 1, new TGeoTranslation(0., 0., 0.));
  
  dz = (kZFilterIn + kZFilterOut) / 2.;
  gMC->Gspos("YMuonFilterM", 1, "YOUT2", 0., 0., - dz, 0, "ONLY");

//
// Outer Pb Cone
//  
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
      par0[11]  = 36.9;  
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
      par0[23]  = 36.9;  // recess erice2000

      par0[24]  = -dz + kZch52;
      par0[25]  = 30.;
      par0[26]  = par0[23];

      par0[27]  = -dz + kZch52;
      par0[28]  = 30.;
      par0[29]  = 30.+(kZch52+4.-kZConeE)*TMath::Tan(kThetaOpenPbO);
// end of cone
      par0[30]  = +dl;
      par0[31]  = 30.;
      par0[32]  = par0[29];
//
      gMC->Gsvolu("YOPB", "PCON", idtmed[kPb+40], par0, 33);
      Float_t dzs = -kzLength + (kZch32-zstart) + dl;
      gMC->Gspos("YOPB", 1, "YMOT", 0., 0., dzs, 0, "ONLY");

//
// Steel envelope
//
      par0[ 0]  = 0.;
      par0[ 1]  = 360.;
      par0[ 2]  = 11.;
  
      par0[ 3]  = -dl;
      par0[ 5]  = 30.+(kZch32-kZConeE)*TMath::Tan(kThetaOpenPbO);
      par0[ 4]  = par0[ 5] - 4.;

//    4th station

      par0[ 6]  = -dz + kZch41 - 4.;
      par0[ 8]  = 30.+(kZch41-4.-kZConeE)*TMath::Tan(kThetaOpenPbO);
      par0[ 7]  = par0[ 8] -4.;

      par0[ 9]  = -dz + kZch41 - 4.;
      par0[11]  = par0[8];  
      par0[10]  = 33.5;

      par0[12]  = -dz + kZch41;
      par0[14]  = 30.+(kZch41-kZConeE)*TMath::Tan(kThetaOpenPbO);  
      par0[13]  = 33.5;

      par0[15]  = -dz + kZch41;
      par0[17]  = 36.9;  
      par0[16]  = 32.9;

//    5th station

      par0[18]  = -dz + kZch51;
      par0[20]  = 36.9;
      par0[19]  = 32.9;

      par0[21]  = -dz + kZch52;
      par0[23]  = 36.9;
      par0[22]  = 32.9;

      par0[24]  = -dz + kZch52;
      par0[26]  = 30.+(kZch52-kZConeE)*TMath::Tan(kThetaOpenPbO);
      par0[25]  = 33.5;

      par0[27]  = -dz + kZch52 + 4.;
      par0[29]  = 30.+(kZch52+4.-kZConeE)*TMath::Tan(kThetaOpenPbO);
      par0[28]  = 33.5;

      par0[30]  = -dz + kZch52 + 4.;
      par0[32]  = 30.+(kZch52+4.-kZConeE)*TMath::Tan(kThetaOpenPbO);
      par0[31]  = par0[32] - 4.;

      par0[33]  = +dl;
      par0[35]  = par0[32];
      par0[34]  = par0[31];

      gMC->Gsvolu("YOSE",    "PCON", idtmed[kSteel], par0, 36);
      gMC->Gspos ("YOSE", 1, "YOPB", 0., 0., 0., 0, "ONLY");
//
//    Extra Tungsten shield close to stations 1 and 2
//
      TGeoRotation* rot000 = new TGeoRotation("rot000",  90.,   0., 90.,  90., 0., 0.);
      TGeoRotation* rot090 = new TGeoRotation("rot090",  90.,  90., 90., 180., 0., 0.);
      TGeoRotation* rot180 = new TGeoRotation("rot180",  90., 180., 90., 270., 0., 0.);
      TGeoRotation* rot270 = new TGeoRotation("rot270",  90., 270., 90.,   0., 0., 0.);
      TGeoVolume* mother = gGeoManager->GetVolume("YMOT");
      TGeoVolumeAssembly* assembly = new TGeoVolumeAssembly("YASS");
      assembly->AddNode(mother, 1, new TGeoTranslation(0., 0., 0.));
      TGeoVolumeAssembly* extraShield1 = new TGeoVolumeAssembly("YCRE");
      TGeoVolumeAssembly* extraShield2 = new TGeoVolumeAssembly("YCRF");

///////////////////////////////////
//                               //
// Recess Station 1              //
//                               //
///////////////////////////////////


///////////////////////////////////
//    FA W-Ring 2                //
//    Drawing ALIP2A__0220       //
///////////////////////////////////
      Float_t faWring2Rinner  = 15.40;
      Float_t faWring2Router  = 18.40;
      Float_t faWring2HWidth  =  3.75;
      Float_t faWring2Cutoffx =  3.35;
      Float_t faWring2Cutoffy =  3.35;
      TGeoTubeSeg* shFaWring2a  = new TGeoTubeSeg(faWring2Rinner, faWring2Router, faWring2HWidth, 0., 90.);
      shFaWring2a->SetName("shFaWring2a");
      TGeoBBox* shFaWring2b  = new TGeoBBox(faWring2Router / 2., faWring2Router / 2., faWring2HWidth);
      shFaWring2b->SetName("shFaWring2b");
      TGeoTranslation* trFaWring2b 
	  = new TGeoTranslation("trFaWring2b", faWring2Router / 2. + faWring2Cutoffx, faWring2Router / 2. + faWring2Cutoffy, 0.);
      trFaWring2b->RegisterYourself();
      TGeoCompositeShape*  shFaWring2 = new TGeoCompositeShape("shFaWring2", "(shFaWring2a)*(shFaWring2b:trFaWring2b)");
      TGeoVolume* voFaWring2    = new TGeoVolume("FA_WRING2", shFaWring2, gGeoManager->GetMedium("SHIL_Ni/W3"));

///////////////////////////////////
//    FA W-Ring 3                //
//    Drawing ALIP2A__0219       //
///////////////////////////////////
      Float_t faWring3Rinner  = 15.40;
      Float_t faWring3Router  = 18.40;
      Float_t faWring3HWidth  =  3.75;
      Float_t faWring3Cutoffx =  3.35;
      Float_t faWring3Cutoffy =  3.35;
      TGeoTubeSeg* shFaWring3a  = new TGeoTubeSeg(faWring3Rinner, faWring3Router, faWring3HWidth, 0., 90.);
      shFaWring3a->SetName("shFaWring3a");
      TGeoBBox* shFaWring3b  = new TGeoBBox(faWring3Router / 2., faWring3Router / 2., faWring3HWidth);
      shFaWring3b->SetName("shFaWring3b");
      TGeoTranslation* trFaWring3b 
	  = new TGeoTranslation("trFaWring3b", faWring3Router / 2. + faWring3Cutoffx, faWring3Router / 2. + faWring3Cutoffy, 0.);
      trFaWring3b->RegisterYourself();
      TGeoCompositeShape*  shFaWring3 = new TGeoCompositeShape("shFaWring3", "(shFaWring3a)*(shFaWring3b:trFaWring3b)");
      TGeoVolume* voFaWring3    = new TGeoVolume("FA_WRING3", shFaWring3, gGeoManager->GetMedium("SHIL_Ni/W3"));

///////////////////////////////////
//    FA W-Ring 5                //
//    Drawing ALIP2A__0221       //
///////////////////////////////////
      Float_t faWring5Rinner = 15.40;
      Float_t faWring5Router = 18.67;
      Float_t faWring5HWidth =  1.08;
      TGeoVolume* voFaWring5    = new TGeoVolume("FA_WRING5", 
						   new TGeoTube(faWring5Rinner, faWring5Router, faWring5HWidth), 
						   gGeoManager->GetMedium("SHIL_Ni/W3"));

//
// Position the rings in the assembly 
//      
// Distance between rings
      Float_t faDWrings = 1.92;
//
      dz = - (4. * faWring2HWidth + 4. * faWring3HWidth + 2. * faWring5HWidth + 2. *  faDWrings) / 2.;
      dz +=  faWring2HWidth;
      extraShield1->AddNode(voFaWring2,    1, new TGeoCombiTrans(0., 0., dz, rot090));
      extraShield1->AddNode(voFaWring2,    2, new TGeoCombiTrans(0., 0., dz, rot270));
      dz +=   faWring2HWidth;      dz +=   faDWrings;
      dz +=   faWring3HWidth;
      extraShield1->AddNode(voFaWring3,    1, new TGeoCombiTrans(0., 0., dz, rot000));
      extraShield1->AddNode(voFaWring3,    2, new TGeoCombiTrans(0., 0., dz, rot180));
      dz +=   faWring3HWidth;   
      dz +=   faWring5HWidth;   
      extraShield1->AddNode(voFaWring5,    1, new TGeoTranslation(0., 0., dz));
      dz +=   faWring5HWidth;   
      dz +=   faWring3HWidth;   
      extraShield1->AddNode(voFaWring3,    3, new TGeoCombiTrans(0., 0., dz, rot090));
      extraShield1->AddNode(voFaWring3,    4, new TGeoCombiTrans(0., 0., dz, rot270));
      dz +=   faWring3HWidth;   
      dz +=   faDWrings;
      dz +=   faWring2HWidth;
      extraShield1->AddNode(voFaWring2,    3, new TGeoCombiTrans(0., 0., dz, rot000));
      extraShield1->AddNode(voFaWring2,    4, new TGeoCombiTrans(0., 0., dz, rot180));
      dz +=   faWring2HWidth;

      // assembly->AddNode(extraShield1, 1, new TGeoTranslation(0., 0., -kzLength + 49.7 + dz));
      Float_t dzKeep = dz;


///////////////////////////////////
//                               //
// Recess Station 2              //
//                               //
///////////////////////////////////

///////////////////////////////////
//    SAA1 W-Ring 1              //
//    Drawing ALIP2A__0217       //
///////////////////////////////////
      Float_t saa1Wring1Width  =  5.85;
      TGeoPcon* shSaa1Wring1  = new TGeoPcon(0., 360., 2);
      shSaa1Wring1->DefineSection(0, 0.00           , 20.30, 23.175);
      shSaa1Wring1->DefineSection(1, saa1Wring1Width, 20.30, 23.400);
      TGeoVolume* voSaa1Wring1  =  new TGeoVolume("SAA1_WRING1", shSaa1Wring1, gGeoManager->GetMedium("SHIL_Ni/W3"));

///////////////////////////////////
//    SAA1 W-Ring 2              //
//    Drawing ALIP2A__0055       //
///////////////////////////////////
      Float_t saa1Wring2Rinner  = 20.30;
      Float_t saa1Wring2Router  = 23.40;
      Float_t saa1Wring2HWidth  =  3.75;
      Float_t saa1Wring2Cutoffx =  4.45;
      Float_t saa1Wring2Cutoffy =  4.45;
      TGeoTubeSeg* shSaa1Wring2a  = new TGeoTubeSeg(saa1Wring2Rinner, saa1Wring2Router, saa1Wring2HWidth, 0., 90.);
      shSaa1Wring2a->SetName("shSaa1Wring2a");
      TGeoBBox* shSaa1Wring2b  = new TGeoBBox(saa1Wring2Router / 2., saa1Wring2Router / 2., saa1Wring2HWidth);
      shSaa1Wring2b->SetName("shSaa1Wring2b");
      TGeoTranslation* trSaa1Wring2b 
	  = new TGeoTranslation("trSaa1Wring2b", saa1Wring2Router / 2. + saa1Wring2Cutoffx, saa1Wring2Router / 2. + saa1Wring2Cutoffy, 0.);
      trSaa1Wring2b->RegisterYourself();
      TGeoCompositeShape*  shSaa1Wring2 = new TGeoCompositeShape("shSaa1Wring2", "(shSaa1Wring2a)*(shSaa1Wring2b:trSaa1Wring2b)");
      TGeoVolume* voSaa1Wring2    = new TGeoVolume("SAA1_WRING2", shSaa1Wring2, gGeoManager->GetMedium("SHIL_Ni/W3"));

///////////////////////////////////
//    SAA1 W-Ring 3              //
//    Drawing ALIP2A__0216       //
///////////////////////////////////

      Float_t saa1Wring3Rinner  = 20.30;
      Float_t saa1Wring3Router  = 23.40;
      Float_t saa1Wring3HWidth  =  3.75;
      Float_t saa1Wring3Cutoffx =  4.50;
      Float_t saa1Wring3Cutoffy =  4.40;
      TGeoTubeSeg* shSaa1Wring3a  = new TGeoTubeSeg(saa1Wring3Rinner, saa1Wring3Router, saa1Wring3HWidth, 0., 90.);
      shSaa1Wring3a->SetName("shSaa1Wring3a");
      TGeoBBox* shSaa1Wring3b  = new TGeoBBox(saa1Wring3Router / 2., saa1Wring3Router / 2., saa1Wring3HWidth);
      shSaa1Wring3b->SetName("shSaa1Wring3b");
      TGeoTranslation* trSaa1Wring3b 
	  = new TGeoTranslation("trSaa1Wring3b", saa1Wring3Router / 2. + saa1Wring3Cutoffx, saa1Wring3Router / 2. + saa1Wring3Cutoffy, 0.);
      trSaa1Wring3b->RegisterYourself();
      TGeoCompositeShape*  shSaa1Wring3 = new TGeoCompositeShape("shSaa1Wring3", "(shSaa1Wring3a)*(shSaa1Wring3b:trSaa1Wring3b)");
      TGeoVolume* voSaa1Wring3    = new TGeoVolume("SAA1_WRING3", shSaa1Wring3, gGeoManager->GetMedium("SHIL_Ni/W3"));

///////////////////////////////////
//    SAA1 W-Ring 4              //
//    Drawing ALIP2A__0215       //
///////////////////////////////////
      Float_t saa1Wring4Width  =  5.85;
      TGeoPcon* shSaa1Wring4  = new TGeoPcon(0., 360., 5);
      shSaa1Wring4->DefineSection(0, 0.00, 20.30, 23.40);
      shSaa1Wring4->DefineSection(1, 1.00, 20.30, 23.40);
      shSaa1Wring4->DefineSection(2, 1.00, 20.30, 24.50);      
      shSaa1Wring4->DefineSection(3, 4.85, 20.30, 24.80);
      shSaa1Wring4->DefineSection(4, 5.85, 24.10, 24.80);
      TGeoVolume* voSaa1Wring4  =  new TGeoVolume("SAA1_WRING4", shSaa1Wring4, gGeoManager->GetMedium("SHIL_Ni/W3"));

///////////////////////////////////
//    SAA1 W-Ring 5              //
//    Drawing ALIP2A__0218       //
///////////////////////////////////
      Float_t saa1Wring5Rinner = 20.30;
      Float_t saa1Wring5Router = 23.40;
      Float_t saa1Wring5HWidth =  0.85;
      TGeoVolume* voSaa1Wring5    = new TGeoVolume("SAA1_WRING5", 
						   new TGeoTube(saa1Wring5Rinner, saa1Wring5Router, saa1Wring5HWidth), 
						   gGeoManager->GetMedium("SHIL_Ni/W3"));
//
// Position the rings in the assembly 
//      
// Distance between rings
      Float_t saa1DWrings = 2.6;
//
      dz = - (saa1Wring1Width + 6. * saa1Wring2HWidth + 2. * saa1Wring3HWidth + saa1Wring4Width + 2. * saa1Wring5HWidth + 2. * saa1DWrings) / 2.;
      extraShield2->AddNode(voSaa1Wring1,    1, new TGeoTranslation(0., 0., dz));
      dz +=   saa1Wring1Width;
      dz +=   saa1Wring2HWidth;   
      extraShield2->AddNode(voSaa1Wring2,    1, new TGeoCombiTrans(0., 0., dz, rot000));
      extraShield2->AddNode(voSaa1Wring2,    2, new TGeoCombiTrans(0., 0., dz, rot180));
      dz +=   saa1Wring2HWidth;   
      dz +=   saa1DWrings;
      dz +=   saa1Wring2HWidth;   
      extraShield2->AddNode(voSaa1Wring2,    3, new TGeoCombiTrans(0., 0., dz, rot090));
      extraShield2->AddNode(voSaa1Wring2,    4, new TGeoCombiTrans(0., 0., dz, rot270));
      dz +=   saa1Wring2HWidth;   
      dz +=   saa1Wring5HWidth;   
      extraShield2->AddNode(voSaa1Wring5,    1, new TGeoTranslation(0., 0., dz));
      dz +=   saa1Wring5HWidth;   
      dz +=   saa1Wring2HWidth;   
      extraShield2->AddNode(voSaa1Wring2,    5, new TGeoCombiTrans(0., 0., dz, rot000));
      extraShield2->AddNode(voSaa1Wring2,    6, new TGeoCombiTrans(0., 0., dz, rot180));
      dz +=   saa1Wring2HWidth;   
      dz +=   saa1DWrings;
      dz +=   saa1Wring3HWidth;   
      extraShield2->AddNode(voSaa1Wring3,    1, new TGeoCombiTrans(0., 0., dz, rot090));
      extraShield2->AddNode(voSaa1Wring3,    2, new TGeoCombiTrans(0., 0., dz, rot270));
      dz +=   saa1Wring3HWidth;   
      extraShield2->AddNode(voSaa1Wring4,    1, new TGeoTranslation(0., 0., dz));
      dz +=   saa1Wring4Width;   
      //assembly->AddNode(extraShield2, 1, new TGeoTranslation(0., 0., -kzLength + (kZch21 - zstart) + dz));
      
      TGeoRotation* rotxz = new TGeoRotation("rotxz",  90., 0., 90., 90., 180., 0.);
      TGeoVolume* yout1 =  gGeoManager->GetVolume("YOUT1");
      yout1->AddNode(extraShield1, 1, new TGeoCombiTrans(0., 0., -zstart - ( 49.7 + dzKeep), rotxz));
      yout1->AddNode(extraShield2, 1, new TGeoCombiTrans(0., 0., -zstart - (kZch21 - zstart + dz ), rotxz));

      TGeoVolume* top =  gGeoManager->GetVolume("ALIC");
      top->AddNode(assembly, 1, new TGeoCombiTrans(0., 0., -zstart - kzLength, rotxz));
  }
}

void AliSHILv2::Init()
{
  //
  // Initialise the muon shield after it has been built
  //
  Int_t i;
  //
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" SHILv2_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    // Here the SHIL initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}
