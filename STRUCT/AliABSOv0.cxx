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
Revision 1.7  2000/10/02 21:28:15  fca
Removal of useless dependecies via forward declarations

Revision 1.6  2000/06/15 09:40:31  morsch
Obsolete typedef keyword removed

Revision 1.5  2000/06/12 19:39:01  morsch
New structure of beam pipe and heating jacket.

Revision 1.4  2000/04/03 08:13:40  fca
Introduce extra scope for non ANSI compliant C++ compilers

Revision 1.3  2000/01/18 17:49:56  morsch
Serious overlap of ABSM with shield corrected
Small error in ARPB parameters corrected

Revision 1.2  2000/01/13 11:23:59  morsch
Last layer of Pb outer angle corrected

Revision 1.1  2000/01/12 15:39:30  morsch
Standard version of ABSO

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
#include "AliMC.h"
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
    
    //
    //

    enum {kC=1605, kAl=1608, kFe=1609, kCu=1610, kW=1611, kPb=1612,
		  kNiCuW=1620, kVacuum=1615, kAir=1614, kConcrete=1616,
		  kPolyCH2=1617, kSteel=1609, kInsulation=1613, kPolyCc=1619};	  
    
    Int_t *idtmed = fIdtmed->GetArray()-1599;
    
    Float_t par[24], cpar[5], cpar0[5], pcpar[12], tpar[3], tpar0[3]; 
    Float_t dz;
#include "ABSOSHILConst.h"
#include "ABSOConst.h"
//
// Structure of Tracking Region
//
  Float_t dzFe = 25.;

// 3 < theta < 9
    fNLayers[0] = 5; 
    fMLayers[0][0]  = kAir;              fZLayers[0][0] = zAbsStart;
    fMLayers[0][1]  = kC;                fZLayers[0][1] = zAbsCc;             
    fMLayers[0][2]  = kConcrete;         fZLayers[0][2] = zRear-dRear-dzFe;
    fMLayers[0][3]  = kFe;               fZLayers[0][3] = zRear-dRear;
    fMLayers[0][4]  = kCu;               fZLayers[0][4] = zRear;
// 2 < theta < 3
    fNLayers[1] = 5; 
    fMLayers[1][0] = fMLayers[0][0];      fZLayers[1][0] = fZLayers[0][0];
    fMLayers[1][1] = fMLayers[0][1];      fZLayers[1][1] = fZLayers[0][1];
    fMLayers[1][2] = fMLayers[0][2];      fZLayers[1][2] = fZLayers[0][2];
    fMLayers[1][3] = fMLayers[0][3];      fZLayers[1][3] = fZLayers[0][3];
    fMLayers[1][4] = kNiCuW;              fZLayers[1][4] = fZLayers[0][4];
//    

    Float_t dTube=0.1;                     // tube thickness
    Float_t dInsu=0.5;                     // insulation thickness
    Float_t dEnve=0.1;                     // protective envelope thickness
    Float_t dFree=0.5;                     // clearance thickness


// Mother volume and outer shielding: Pb
  par[0]  = 0.;
  par[1]  = 360.;
  par[2]  = 7.;
    
  par[3]  = -(zRear-zAbsStart)/2.;
  par[4]  = rAbs;
  par[5]  = zAbsStart * TMath::Tan(theta1);

  par[6]  = par[3]+(zNose-zAbsStart);
  par[7]  = rAbs;
  par[8]  = zNose * TMath::Tan(theta1);

  par[9]  = par[3]+(zConeTPC-zAbsStart);
  par[10] = rAbs;
  par[11] = par[8] + (par[9] - par[6]) * TMath::Tan(theta2);

  par[12]  = par[3]+(zOpen-zAbsStart);
  par[13] = rAbs;
  par[14] = par[11] + (par[12] - par[9]) * TMath::Tan(accMax);

  par[15] = par[3]+(zRear-dRear-zAbsStart);
  par[16] = rAbs   + (par[15] - par[12]) * TMath::Tan(thetaOpen1) ;
  par[17] = par[14] + (par[15] - par[12]) * TMath::Tan(accMax);

  par[18] = par[3]+(zRear-dRear-zAbsStart);
  par[19] = (zRear-dRear) * TMath::Tan(accMin);
  par[20] = par[14] + (par[18] - par[12]) * TMath::Tan(accMax);

  par[21] = -par[3];
  par[22] =  zRear* TMath::Tan(accMin);
  par[23] = par[20] + (par[21] - par[18]) * TMath::Tan(accMax);
  gMC->Gsvolu("ABSS", "PCON", idtmed[kPb], par, 24);
  { // Begin local scope for i
      for (Int_t i=4; i<18; i+=3) par[i]  = 0;
  } // End local scope for i
  gMC->Gsvolu("ABSM", "PCON", idtmed[kVacuum+40], par, 24);
  gMC->Gspos("ABSS", 1, "ABSM", 0., 0., 0., 0, "ONLY");

//
// Steel envelope
//
  par[4] = par[5] -dSteel;
  par[7] = par[8] -dSteel;
  par[10]= par[11]-dSteel;  
  par[13]= par[14]-dSteel;  
  par[16]= par[17]-dSteel;  
  par[19]= par[20]-dSteel;  
  par[22]= par[23]-dSteel;  
  gMC->Gsvolu("ABST", "PCON", idtmed[kSteel], par, 24);
  gMC->Gspos("ABST", 1, "ABSS", 0., 0., 0., 0, "ONLY");
//
// Polyethylene shield
// 
  cpar[0] = (zRear - zConeTPC) / 2.;
  cpar[1] = zConeTPC * TMath::Tan(accMax);
  cpar[2] = cpar[1] + dPoly;
  cpar[3] = zRear * TMath::Tan(accMax);
  cpar[4] = cpar[3] + dPoly;
  gMC->Gsvolu("APOL", "CONE", idtmed[kPolyCH2+40], cpar, 5);
  dz = (zRear-zAbsStart)/2.-cpar[0];
  gMC->Gspos("APOL", 1, "ABSS", 0., 0., dz, 0, "ONLY");

//
// Tungsten nose to protect TPC
// 
  cpar[0] = (zNose - zAbsStart) / 2.;
  cpar[1] = zAbsStart * TMath::Tan(accMax);
  cpar[2] = zAbsStart * TMath::Tan(theta1)-dSteel;
  cpar[3] = zNose * TMath::Tan(accMax);
  cpar[4] = zNose * TMath::Tan(theta1)-dSteel;
  gMC->Gsvolu("ANOS", "CONE", idtmed[kW], cpar, 5);
//
  dz = -(zRear-zAbsStart)/2.+cpar[0];
  gMC->Gspos("ANOS", 1, "ABSS", 0., 0., dz, 0, "ONLY");
//
// Tungsten inner shield
//
  Float_t zW=zTwoDeg+.1;
  Float_t dZ = zW+(zRear-dRear-zW)/2.;
  //
  pcpar[0]  = 0.;
  pcpar[1]  = 360.;
  pcpar[2]  = 3.;
  pcpar[3]  = zW-dZ;
  pcpar[4]  = rAbs;
  pcpar[5]  = zW * TMath::Tan(accMin);
  pcpar[6]  = zOpen-dZ;
  pcpar[7]  = rAbs;
  pcpar[8]  = zOpen * TMath::Tan(accMin);
  pcpar[9]  = zRear-dRear-dZ;
  pcpar[10] = rAbs+(zRear-dRear-zOpen) * TMath::Tan(thetaOpen1);
  pcpar[11] = (zRear-dRear) * TMath::Tan(accMin);
  
  gMC->Gsvolu("AWIN", "PCON", idtmed[kNiCuW+40], pcpar, 12);
  //
  dz=(zW+zRear-dRear)/2-(zAbsStart+zRear)/2.;
  gMC->Gspos("AWIN", 1, "ABSS", 0., 0., dz, 0, "ONLY");

  //     Inner tracking region
  //
  //     mother volume: Cu
  //
  pcpar[0]  = 0.;
  pcpar[1]  = 360.;
  pcpar[2]  = 3.;
  pcpar[3]  = -(zRear-zAbsStart)/2.;
  pcpar[4]  = rAbs;
  pcpar[5]  = zAbsStart * TMath::Tan(accMax);
  pcpar[6]  = pcpar[3]+(zTwoDeg-zAbsStart);
  pcpar[7]  = rAbs;
  pcpar[8]  = zTwoDeg * TMath::Tan(accMax);
  pcpar[9]  = -pcpar[3];
  pcpar[10] = zRear * TMath::Tan(accMin);
  pcpar[11] = zRear * TMath::Tan(accMax);
  gMC->Gsvolu("AITR", "PCON", idtmed[fMLayers[0][4]], pcpar, 12);
  //
  // special Pb medium for last 5 cm of Pb
  Float_t zr=zRear-2.-0.001;
  cpar[0] = 1.0;
  cpar[1] = zr * TMath::Tan(thetaR);
  cpar[2] = zr * TMath::Tan(accMax);
  cpar[3] = cpar[1] + TMath::Tan(thetaR) * 2;
  cpar[4] = cpar[2] + TMath::Tan(accMax) * 2;
  gMC->Gsvolu("ARPB", "CONE", idtmed[fMLayers[0][4]], cpar, 5);
  dz=(zRear-zAbsStart)/2.-cpar[0]-0.001;
  gMC->Gspos("ARPB", 1, "AITR", 0., 0., dz, 0, "ONLY");
  //
  //     concrete cone: concrete 
  //
  pcpar[9]  = pcpar[3]+(zRear-dRear-zAbsStart);
  pcpar[10] = (zRear-dRear) * TMath::Tan(accMin);
  pcpar[11] = (zRear-dRear) * TMath::Tan(accMax);
  gMC->Gsvolu("ACON", "PCON", idtmed[fMLayers[0][2]+40], pcpar, 12);
  gMC->Gspos("ACON", 1, "AITR", 0., 0., 0., 0, "ONLY");
//
//    Fe Cone 
//
  zr = zRear-dRear-dzFe;
  cpar[0]  = dzFe/2.;
  cpar[1] = zr * TMath::Tan(accMin);
  cpar[2] = zr * TMath::Tan(accMax);
  cpar[3] = cpar[1] + TMath::Tan(thetaR) * dzFe;
  cpar[4] = cpar[2] + TMath::Tan(accMax) * dzFe;
  gMC->Gsvolu("ACFE", "CONE",idtmed[fMLayers[0][3]], cpar, 5);

  dz = (zRear-zAbsStart)/2.-dRear-dzFe/2.;

  gMC->Gspos("ACFE", 1, "ACON", 0., 0., dz, 0, "ONLY");

  
  //
  //
  //     carbon cone: carbon
  //
  pcpar[9]  = pcpar[3]+(zAbsCc-zAbsStart);
  pcpar[10]  = zAbsCc * TMath::Tan(accMin);
  pcpar[11]  = zAbsCc * TMath::Tan(accMax);
  gMC->Gsvolu("ACAR", "PCON", idtmed[fMLayers[0][1]+40], pcpar, 12);
  gMC->Gspos("ACAR", 1, "ACON", 0., 0., 0., 0, "ONLY");
 //
 //     carbon cone outer region
 //
  cpar[0]  = 10.;
  cpar[1]  = rAbs;
  cpar[2]  = zAbsStart* TMath::Tan(accMax);
  cpar[3]  = rAbs;
  cpar[4]  = cpar[2]+2. * cpar[0] * TMath::Tan(accMax);

  gMC->Gsvolu("ACAO", "CONE", idtmed[fMLayers[0][1]], cpar, 5);
  dz=-(zRear-zAbsStart)/2.+cpar[0];
  gMC->Gspos("ACAO", 1, "ACAR", 0., 0., dz, 0, "ONLY");
  //
  //     inner W shield
  Float_t epsi=0.;
  Float_t repsi=1.;
  
  zr=zRear-(dRear-epsi);
  cpar[0] = (dRear-epsi)/2.;
  cpar[1] = zr * TMath::Tan(accMin);
  cpar[2] = zr * TMath::Tan(thetaR*repsi);
  cpar[3] = cpar[1] + TMath::Tan(accMin) * (dRear-epsi);
  cpar[4] = cpar[2] + TMath::Tan(thetaR*repsi) * (dRear-epsi);
  gMC->Gsvolu("ARW0", "CONE", idtmed[fMLayers[1][4]+40], cpar, 5);
  dz=(zRear-zAbsStart)/2.-cpar[0];
  gMC->Gspos("ARW0", 1, "AITR", 0., 0., dz, 0, "ONLY");
  //
  // special W medium for last 5 cm of W
  zr=zRear-5;
  cpar[0] = 2.5;
  cpar[1] = zr * TMath::Tan(accMin);
  cpar[2] = zr * TMath::Tan(thetaR*repsi);
  cpar[3] = cpar[1] + TMath::Tan(accMin) * 5.;
  cpar[4] = cpar[2] + TMath::Tan(thetaR*repsi) * 5.;
  gMC->Gsvolu("ARW1", "CONE", idtmed[fMLayers[1][4]+20], cpar, 5);
  dz=(dRear-epsi)/2.-cpar[0];
  gMC->Gspos("ARW1", 1, "ARW0", 0., 0., dz, 0, "ONLY");
  //
  // Cu
  Float_t drMin=TMath::Tan(thetaR) * 5;
  Float_t drMax=TMath::Tan(accMax) * 5;
  gMC->Gsvolu("ARPE", "CONE", idtmed[fMLayers[0][4]], cpar, 0);
  cpar[0]=2.5;
  { // Begin local scope for i
      for (Int_t i=0; i<3; i++) {
	  zr=zRear-dRear+5+i*10.;
	  cpar[1] = zr * TMath::Tan(thetaR);
	  cpar[2] = zr * TMath::Tan(accMax);
	  cpar[3] = cpar[1] + drMin;
	  cpar[4] = cpar[2] + drMax;
	  dz=(zRear-zAbsStart)/2.-cpar[0]-5.-(2-i)*10;
	  gMC->Gsposp("ARPE", i+1, "AITR", 0., 0., dz, 0, "ONLY",cpar,5);
      }
  } // End local scope for i
  gMC->Gspos("AITR", 1, "ABSS", 0., 0., 0., 0, "ONLY");	
  dz = (zRear-zAbsStart)/2.+zAbsStart;
  gMC->Gspos("ABSM", 1, "ALIC", 0., 0., dz, 0, "ONLY");	
//
//
// vacuum system
//
// pipe and heating jackets
//
//
// cylindrical piece
  tpar0[2]=(zOpen-zAbsStart)/2;
  tpar0[0]=rVacu;
  tpar0[1]=rAbs;
  gMC->Gsvolu("AV11", "TUBE", idtmed[kSteel+40], tpar0, 3);
//
// insulation

  tpar[2]=tpar0[2];
  tpar[0]=rVacu+dTube;
  tpar[1]=tpar[0]+dInsu;
  gMC->Gsvolu("AI11", "TUBE", idtmed[kInsulation+40], tpar, 3);
  gMC->Gspos("AI11", 1, "AV11", 0., 0., 0., 0, "ONLY"); 
//
// clearance 
  tpar[0]=tpar[1]+dEnve;
  tpar[1]=tpar[0]+dFree;
  gMC->Gsvolu("AP11", "TUBE", idtmed[kAir+40], tpar, 3);
  gMC->Gspos("AP11", 1, "AV11", 0., 0., 0., 0, "ONLY"); 
//
  dz=-(zRear-zAbsStart)/2.+tpar0[2];
  gMC->Gspos("AV11", 1, "ABSM", 0., 0., dz, 0, "ONLY"); 
//
// conical piece

  cpar0[0]=(zRear-dRear-zOpen)/2;
  cpar0[1]=rVacu-0.05;
  cpar0[2]=rAbs;
  Float_t dR=2.*cpar0[0]*TMath::Tan(thetaOpen1);
  cpar0[3]=cpar0[1]+dR;
  cpar0[4]=cpar0[2]+dR;
  gMC->Gsvolu("AV21", "CONE", idtmed[kSteel+40], cpar0, 5);
  dTube+=0.05;

//
// insulation
  cpar[0]=cpar0[0];
  cpar[1]=cpar0[1]+dTube;
  cpar[2]=cpar0[1]+dTube+dInsu;
  cpar[3]=cpar0[3]+dTube;
  cpar[4]=cpar0[3]+dTube+dInsu;
  gMC->Gsvolu("AI21", "CONE", idtmed[kInsulation+40], cpar, 5);
  gMC->Gspos("AI21", 1, "AV21", 0., 0., 0., 0, "ONLY"); 
//
// clearance
  cpar[1]=cpar0[1]+dTube+dInsu+dEnve;
  cpar[2]=rAbs;
  cpar[3]=cpar0[1]+dTube+dInsu+dEnve+dR;
  cpar[4]=rAbs+dR;

  gMC->Gsvolu("AP21", "CONE", idtmed[kAir+40], cpar, 5);
  gMC->Gspos("AP21", 1, "AV21", 0., 0., 0., 0, "ONLY"); 
  
  dz=(zRear-zAbsStart)/2.-cpar0[0]-dRear;
  gMC->Gspos("AV21", 1, "ABSM", 0., 0., dz, 0, "ONLY"); 
//
// Support cone 

  par[0]  = 0.;
  par[1]  = 360.;
  par[2]  = 4.;
    
  par[3]  = zRear;
  par[4]  = 100.;
  par[5]  = 170.;
  
  par[6]  = zRear+2.;
  par[7]  = 100.;
  par[8]  = 170.;

  par[9]  = zRear+2.;
  par[10] = 168.;
  par[11] = 170.;

  par[12]  = 600.;
  par[13] = 168.;
  par[14] = 170.;
  

  gMC->Gsvolu("ASSS", "PCON", idtmed[kSteel], par, 25);
  gMC->Gspos("ASSS", 1, "ALIC", 0., 0., 0., 0, "ONLY");


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
 









