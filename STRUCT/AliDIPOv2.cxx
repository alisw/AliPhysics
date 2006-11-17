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
//  Magnetic Dipole version 1                                                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliDIPOv2Class.gif">
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

#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliDIPOv2.h"
#include "AliMagF.h"
#include "AliRun.h"
 
ClassImp(AliDIPOv2)
 
//_____________________________________________________________________________
AliDIPOv2::AliDIPOv2() 
{
  //
  // Last design of magnetic dipole version 2
  //
}
 
//_____________________________________________________________________________
AliDIPOv2::AliDIPOv2(const char *name, const char *title)
  : AliDIPO(name,title)
{
  //
  // Standard constructor for the magnetic dipole version 2    
}

void AliDIPOv2::CreateGeometry()
{
  //
  // Creation of the geometry of the magnetic DIPOLE version 2
  //

    CreateSpectrometerDipole();
    CreateCompensatorDipole();
}


//_____________________________________________________________________________
void AliDIPOv2::CreateSpectrometerDipole()
{
  //
  // Creation of the geometry of the magnetic DIPOLE version 2
  //

  Float_t cpar[5], tpar[18], ypar[12];
  Float_t dz, dx, dy;
  Int_t idrotm[1899];
  Float_t accMax, the1, phi1, the2, phi2, the3, phi3;
  
  Int_t *idtmed = fIdtmed->GetArray()-1799;

//  const Int_t kCoil = 1813;
//  const Int_t kCable= 1811;

  const Int_t kCoil = 1808;
  const Int_t kCable= 1808;
  
  accMax = 9.;   // ANGLE POLAIRE MAXIMUM 

  //       DIPOLE MAGNET 
  const Float_t kZDipole = 975; 

  tpar[ 0] = 0.; 
  tpar[ 1] = 360.;
  tpar[ 2] = 4.; 
  //
  tpar[ 3] = -250.55 + kZDipole;
  tpar[ 4] =  30.1;
  tpar[ 5] = 527.34; 
  //
  tpar[ 6] = 37. + kZDipole;
  tpar[ 7] =  30.1;
  tpar[ 8] = 527.34;
  //
  tpar[ 9] = 37. + kZDipole;
  tpar[10] = tpar[ 9] * TMath::Tan(2. * TMath::Pi() / 180.);
  tpar[11] = 527.34;
  //
  tpar[12] = 260.55 + kZDipole;
  tpar[13] = tpar[12] * TMath::Tan(2. * TMath::Pi() / 180.);
  tpar[14] = 527.34;
  gMC->Gsvolu("DDIP", "PCON", idtmed[1874], tpar, 15);
  //
  //       Coils 
  // air - m.f. 
  cpar[0] = 207.;
  cpar[1] = 274.;
  cpar[2] = 37.65;
  cpar[3] = 119.;
  cpar[4] = 241.; 
  //   coil - high cuts
  gMC->Gsvolu("DC1 ", "TUBS", idtmed[kCoil+40], cpar, 5);
  cpar[3] = -61.;
  cpar[4] = 61.;
  gMC->Gsvolu("DC2 ", "TUBS", idtmed[kCoil+40], cpar, 5);

  //  coil - low cuts cuts
  cpar[0] = 207.;
//  cpar[1] = cpar[0] + 10.;
  cpar[1] = 217;
  cpar[3] = 119.;
  cpar[4] = 241.;

  gMC->Gsvolu("DC3 ", "TUBS", idtmed[kCoil], cpar, 5);
  cpar[0] = 207.; 
  cpar[1] = 217; 
  cpar[3] = -61.;
  cpar[4] = 61.;
  gMC->Gsvolu("DC4 ", "TUBS", idtmed[kCoil], cpar, 5);

  gMC->Gspos("DC3 ", 1, "DC1 ", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("DC4 ", 1, "DC2 ", 0., 0., 0., 0, "ONLY");

//  dz =  37.65 - 243.55
  dz = -205.9-2.45;
  dx = 5.;
  gMC->Gspos("DC1 ", 1, "DDIP",  dx, 0.,  dz+kZDipole, 0, "ONLY");
  gMC->Gspos("DC1 ", 2, "DDIP",  dx, 0., -dz+kZDipole, 0, "ONLY");
  gMC->Gspos("DC2 ", 1, "DDIP", -dx, 0.,  dz+kZDipole, 0, "ONLY");
  gMC->Gspos("DC2 ", 2, "DDIP", -dx, 0., -dz+kZDipole, 0, "ONLY");
  the1 = 180.;
  phi1 = 0.;
  the2 = 90.;
  phi2 = 151.;
  the3 = 90.;
  phi3 = 61.;
  AliMatrix(idrotm[1800], the1, phi1, the2, phi2, the3, phi3);
  phi2 = 29.;  //90-61
  the3 = -90.;
  phi3 = -61.;
  AliMatrix(idrotm[1801], the1, phi1, the2, phi2, the3, phi3);
  the1 = 0.;
  phi1 = 0.;
  the2 = 90.;
  phi2 = 151.;
  the3 = 90.;
  phi3 = 61.;
  AliMatrix(idrotm[1802], the1, phi1, the2, phi2, the3, phi3);
  phi2 = 29.;
  the3 = -90.;
  phi3 = -61.;
  AliMatrix(idrotm[1803], the1, phi1, the2, phi2, the3, phi3);

  cpar[0] = 25.;
  cpar[1] = 100.3; //25+75.3
  cpar[2] = 33.5;
  cpar[3] = 270.;
  cpar[4] = 360.;
//*  coil high cuts
  gMC->Gsvolu("DC11", "TUBS", idtmed[kCoil+40], cpar, 5);

  dx = TMath::Sin(30.5*kDegrad) * -(207.+33.5)+5./TMath::Sin(30.5*kDegrad); 
  dy = TMath::Cos(30.5*kDegrad) * -(207.+33.5);  
  dz = cpar[1] - 243.55-2.45;
  gMC->Gspos("DC11", 1, "DDIP",  dx, dy,  dz+kZDipole, idrotm[1800], "ONLY");
  gMC->Gspos("DC11", 2, "DDIP",  dx, dy, -dz+kZDipole, idrotm[1802], "ONLY");
  gMC->Gspos("DC11", 3, "DDIP", -dx, dy,  dz+kZDipole, idrotm[1801], "ONLY");
  gMC->Gspos("DC11", 4, "DDIP", -dx, dy, -dz+kZDipole, idrotm[1803], "ONLY");



//* ... higher cuts
  cpar[0] = 25.;
  cpar[1] = 100.3; //25+75.3
  cpar[2] = 33.5;
  cpar[3] = 0.;
  cpar[4] = 90.;
//*  coil high cuts
  gMC->Gsvolu("DC12", "TUBS", idtmed[kCoil+40], cpar, 5);

  dx = TMath::Sin(30.5*kDegrad) * -(207.+33.5)+5./TMath::Sin(30.5*kDegrad); 
  dy = TMath::Cos(30.5*kDegrad) *(207.+33.5);  
  dz = cpar[1] - 243.55-2.45;
  gMC->Gspos("DC12", 1, "DDIP",  dx, dy,  dz+kZDipole, idrotm[1801], "ONLY");
  gMC->Gspos("DC12", 2, "DDIP",  dx, dy, -dz+kZDipole, idrotm[1803], "ONLY");
  gMC->Gspos("DC12", 3, "DDIP", -dx, dy,  dz+kZDipole, idrotm[1800], "ONLY");
  gMC->Gspos("DC12", 4, "DDIP", -dx, dy, -dz+kZDipole, idrotm[1802], "ONLY");

  the1 = 90.;
  phi1 = 61.;
  the2 = 90.;
  phi2 = 151.;
  the3 = 0.;
  phi3 = 0.;
  AliMatrix(idrotm[1804], the1, phi1, the2, phi2, the3, phi3);
  the1 = 90.;
  phi1 = -61.;
  the2 = 90.;
  phi2 = -151.;
  AliMatrix(idrotm[1805], the1, phi1, the2, phi2, the3, phi3);
  the1 = 90.;
  phi1 = 119.; //180 -61
  the2 = 90.;
  phi2 = 209.; //270-61
  AliMatrix(idrotm[1806], the1, phi1, the2, phi2, the3, phi3);
  the1 = 90.;
  phi1 = -119.;
  the2 = 90.;
  phi2 = -209.;
  AliMatrix(idrotm[1807], the1, phi1, the2, phi2, the3, phi3); 

//*  coil - high cuts

  tpar[0] = 37.65;
  tpar[1] = 33.5;
  tpar[2] = 145.5;
  gMC->Gsvolu("DL1 ", "BOX ", idtmed[kCoil+40], tpar, 3);

// coil - low cuts

  tpar[0] = 5.;
  dx = 37.65  - 5.;  
  gMC->Gsvolu("DL2 ", "BOX ", idtmed[kCoil], tpar, 3);
  gMC->Gspos("DL2 ", 1, "DL1 ", dx, 0., 0., 0, "ONLY");

  dx =-53.62;
  dy =-241.26819;
  dz = 0.; 
  gMC->Gspos("DL1 ", 1, "DDIP", dx,  dy, dz+kZDipole, idrotm[1804], "ONLY");
  gMC->Gspos("DL1 ", 2, "DDIP", dx, -dy, dz+kZDipole, idrotm[1805], "ONLY");
  gMC->Gspos("DL1 ", 3, "DDIP",-dx,  dy, dz+kZDipole, idrotm[1806], "ONLY");
  gMC->Gspos("DL1 ", 4, "DDIP",-dx, -dy, dz+kZDipole, idrotm[1807], "ONLY");

  // Contactor

 //  high cuts

  //Steel outer face planes

  cpar[0] = 207.-18.6;
  cpar[1] = 274.+18.6;
  cpar[2] = 1.;
  cpar[3] = -50.;
  cpar[4] = 50.; 
 
  gMC->Gsvolu("DCO1", "TUBS", idtmed[1818], cpar, 5);
  dx = -5.;
  dz = 168.25-1.5-1.;
  gMC->Gspos("DCO1", 1, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");
  dz = 243.55+4.5+1.5+1.;
  gMC->Gspos("DCO1", 2, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");
  
  // 9.06.2000

  //  cpar[0] = 207.-18.6;
  //  cpar[1] = 274.+18.6;
  // cpar[2] = 1.;
  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 
 
  gMC->Gsvolu("DCO2", "TUBS", idtmed[1818], cpar, 5);
  dx = +5.;
  dz = 168.25-1.5-1.;
  gMC->Gspos("DCO2", 1, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");
  dz = 243.55+4.5+1.5+1.;
  gMC->Gspos("DCO2", 2, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");
 


  // Resin face planes

  cpar[0] = 207.;
  cpar[1] = 274.;
  cpar[2] = .75;
  cpar[3] = -50.;
  cpar[4] = 50.; 
 
  gMC->Gsvolu("DCO3", "TUBS", idtmed[1812], cpar, 5);
  dx = -5;
  dz = 168.25-0.75;
  gMC->Gspos("DCO3", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");
  dz = 243.55+4.5+0.75;
  gMC->Gspos("DCO3", 2, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000

  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 
  gMC->Gsvolu("DCO4", "TUBS", idtmed[1812], cpar, 5);
  dx = +5;
  dz = 168.25-0.75;
  gMC->Gspos("DCO4", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");
  dz = 243.55+4.5+0.75 ;
  gMC->Gspos("DCO4", 2, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

 
  // G10 face plane

  cpar[0] = 207.;
  cpar[1] = 274.;
  cpar[2] = 2.25;
  cpar[3] = -50.;
  cpar[4] = 50.; 
 
  gMC->Gsvolu("DCO5", "TUBS", idtmed[1810], cpar, 5);

  dx = -5;
  dz = 243.55+2.25;
  gMC->Gspos("DCO5", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000

  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 

  gMC->Gsvolu("DCO6", "TUBS", idtmed[1810], cpar, 5);

  dx = +5;
  dz = 243.55+2.25;
  gMC->Gspos("DCO6", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  //Steel supported planes

  cpar[0] = 274.+1.5+2.;
  cpar[1] = 274.+18.6;
  cpar[2] = 1.;
  cpar[3] = -50.;
  cpar[4] = 50.;  
 
  gMC->Gsvolu("DCO7", "TUBS", idtmed[1818], cpar, 5);

  dx = -5;
  dz = 168.25+1.;
  gMC->Gspos("DCO7", 1, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000
  cpar[0] = 274.+1.5+2.;
  cpar[1] = 274.+18.6;
  cpar[2] = 1.;
  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 

 
  gMC->Gsvolu("DCO8", "TUBS", idtmed[1818], cpar, 5);

  dx = +5;
  dz = 168.25+1.;
  gMC->Gspos("DCO8", 1, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");

  //

  cpar[0] = 207.- 18.6;
  cpar[1] = 207.- 2.- 1.5;
  cpar[2] = 1.;
  cpar[3] = -50.;
  cpar[4] = 50.; 

  gMC->Gsvolu("DCO9", "TUBS", idtmed[1818], cpar, 5);

  dx = -5;
  dz = 168.25+1.;
  gMC->Gspos("DCO9", 1, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000

  cpar[0] = 207.- 18.6;
  cpar[1] = 207.- 2.- 1.5;
  cpar[2] = 1.;
  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 

  gMC->Gsvolu("DCOA", "TUBS", idtmed[1818], cpar, 5);

  dx = +5;
  dz = 168.25+1.;
  gMC->Gspos("DCOA", 1, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");


  // Sides steel planes

  cpar[0] = 207. - 1.5 -2.;
  cpar[1] = 207. - 1.5;
  cpar[2] = ((243.55+4.5+1.5)-168.25)/2;
  cpar[3] = -50.;
  cpar[4] = 50.; 
 
  gMC->Gsvolu("DCOB", "TUBS", idtmed[1818], cpar, 5);

  cpar[0] = 274. + 1.5;
  cpar[1] = 274. + 1.5 +2.;

  gMC->Gsvolu("DCOC", "TUBS", idtmed[1818], cpar, 5);

  dx=-5.;
  dz = ((243.55+4.5+1.5)+168.25)/2;
  gMC->Gspos("DCOB", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");
  gMC->Gspos("DCOC", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000

  cpar[0] = 207. - 1.5 -2.;
  cpar[1] = 207. - 1.5;
  cpar[2] = ((243.55+4.5+1.5)-168.25)/2;
  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 

  gMC->Gsvolu("DCOD", "TUBS", idtmed[1818], cpar, 5);

  cpar[0] = 274. + 1.5;
  cpar[1] = 274. + 1.5 +2.;

  gMC->Gsvolu("DCOE", "TUBS", idtmed[1818], cpar, 5);

  dx=+5.;
  dz = ((243.55+4.5+1.5)+168.25)/2;
  gMC->Gspos("DCOD", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");
  gMC->Gspos("DCOE", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");


  // Top and bottom resin  planes

  cpar[0] = 207. - 1.5;
  cpar[1] = 207.;
  cpar[2] = ((243.55+4.5+1.5)-168.25)/2;
  cpar[3] = -50.;
  cpar[4] = 50.; 
 
  gMC->Gsvolu("DCOF", "TUBS", idtmed[1812], cpar, 5);

  cpar[0] = 274.;
  cpar[1] = 274. + 1.5;

  gMC->Gsvolu("DCOG", "TUBS", idtmed[1812], cpar, 5);


  dx=-5.;
  dz = ((243.55+4.5+1.5)+168.25)/2;
  gMC->Gspos("DCOF", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");
  gMC->Gspos("DCOG", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000
  cpar[0] = 207. - 1.5;
  cpar[1] = 207.;
  cpar[2] = ((243.55+4.5+1.5)-168.25)/2;

  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 

  gMC->Gsvolu("DCOH", "TUBS", idtmed[1812], cpar, 5);

  cpar[0] = 274.;
  cpar[1] = 274. + 1.5;

  gMC->Gsvolu("DCOI", "TUBS", idtmed[1812], cpar, 5);


  dx=+5.;
  dz = ((243.55+4.5+1.5)+168.25)/2;
  gMC->Gspos("DCOH", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");
  gMC->Gspos("DCOI", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");


  // Aluminum cabels

  cpar[0] = 274. + 1.5  +2.;
  cpar[1] = 274. + 1.5  +2. + 80.;
  cpar[2] = 5.05/2;
  cpar[3] = -24.;
  cpar[4] = 24.; 
 
  gMC->Gsvolu("DCOJ", "TUBS", idtmed[kCable], cpar, 5);

  //  dx = 274. + 1.5  +2. +40.;
  //  dx = 5. + 1.5 +2. +40.;
  //  dx = 5. + 1.5 +2.;
  dx=-5.;
  dz = 168.25 + 5.05 + 5.05/2;
  gMC->Gspos("DCOJ", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  dz = 243.55 - 5.05/2;
  gMC->Gspos("DCOJ", 2, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000

  cpar[3] = 180.-24.;
  cpar[4] = 180.+24.; 

   gMC->Gsvolu("DCOK", "TUBS", idtmed[kCable], cpar, 5);

  //  dx = 274. + 1.5  +2. +40.;
  //  dx = 5. + 1.5 +2. +40.;
  //  dx = 5. + 1.5 +2.;
  dx=+5.;
  dz = 168.25 + 5.05 + 5.05/2;
  gMC->Gspos("DCOK", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  dz = 243.55 - 5.05/2;
  gMC->Gspos("DCOK", 2, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

 
  //   YOKE 

// Top and bottom blocks
  ypar[0] = 298.1; 
  ypar[1] = 69.5;
  ypar[2] = 155.75;

// iron- high cuts
  gMC->Gsvolu("DY1 ", "BOX ", idtmed[1858], ypar, 3);
  ypar[0] = 144.+10.; 
  ypar[1] = 193.3+10.;
  ypar[2] = 5.;
  ypar[3] = 155.75;
  dy = -69.5 + 5.;
// iron- low cuts
  gMC->Gsvolu("DY11", "TRD1", idtmed[1818], ypar, 4);
  gMC->Gspos("DY11", 1, "DY1 ", 0.,  dy, 0., 0, "ONLY");

  dy = 365.5;
  dz = 4.95;
  gMC->Gspos("DY1 ", 1, "DDIP", 0.,  dy, -dz+kZDipole, 0, "ONLY");

  the1 = 270.;
  phi1 = 0.;
  the2 = 270.;
  phi2 = 90.;
  the3 = 0.;
  phi3 = 0.;
  AliMatrix(idrotm[1808], the1, phi1, the2, phi2, the3, phi3);
  gMC->Gspos("DY1 ", 2, "DDIP", 0., -dy, -dz+kZDipole, idrotm[1808] , "ONLY");

// side walls
  //  ypar[0] = 579./2.; 
  ypar[0] = 296.; 
  ypar[1] = 0.;
  ypar[2] = 0.;
  ypar[3] = 155.75;
  ypar[4] = 47.9;
  ypar[5] = 72.55;
  ypar[6] = 4.3058039629;
  // z+ 
  ypar[7] = 155.75;
  ypar[8] = 47.9;
  ypar[9] = 72.55;
  ypar[10] = 4.3058039629;

// iron - high cuts

  gMC->Gsvolu("DY2 ", "TRAP", idtmed[1858], ypar,11);

  ypar[4] = 47.9 -5.;
  ypar[5] = 72.55 -5.;

  ypar[8] = 47.9 -5.;
  ypar[9] = 72.55 -5.;


// iron - low cuts

  gMC->Gsvolu("DY22", "TRAP", idtmed[1818], ypar,11);

  dy = 0.;
  dx = -5.;

  gMC->Gspos("DY22", 1, "DY2 ", dx,  dy, 0., 0, "ONLY");

  the1 = 90.;
  phi1 = 180.;
  the2 = 180.;
  phi2 = 180.;
  the3 = 90.;
  phi3 = 90.;
  AliMatrix(idrotm[1809], the1, phi1, the2, phi2, the3, phi3);

  the1 = 90.;
  phi1 = 0.;
  the2 = 180.;
  phi2 = 0.;
  the3 = 90.;
  phi3 = 90.;
  AliMatrix(idrotm[1810], the1, phi1, the2, phi2, the3, phi3);

  dx = 228.875;
  dz = - 4.95;
  
  gMC->Gspos("DY2 ", 1, "DDIP",  dx, 0.0,  dz+kZDipole, idrotm[1809], "ONLY");
  gMC->Gspos("DY2 ", 2, "DDIP", -dx, 0.0,  dz+kZDipole, idrotm[1810], "ONLY");
  
  AliMatrix(idrotm[1811], 270., 0., 90., 90., 180., 0.);
  gMC->Gspos("DDIP", 1, "ALIC", 0., 0., 0., idrotm[1811], "ONLY");
  gMC->Gsatt("DDIP", "SEEN", 0);
}


void AliDIPOv2::CreateCompensatorDipole()
{
    //
    //  Geometry of the compensator Dipole MBWMD (was MCB @ SPS)
    // 
    Int_t *idtmed = fIdtmed->GetArray()-1799;
    Int_t idrotm[1899];
//
    Float_t pbox[3] = {63., 63., 170.};
    
//  Mother volumes
    gMC->Gsvolu("DCM0", "BOX", idtmed[1814], pbox, 3);
    gMC->Gsatt("DCM0", "SEEN", 0);
    
//
//  Mother volume containing lower coil
    pbox[0] = 58.5/2.;
    pbox[1] = 30.0;
    pbox[2] = 152.5;
    
    gMC->Gsvolu("DCML", "BOX", idtmed[1809], pbox, 3);
//
// Base
    pbox[0] = 62.5;
    pbox[1] = 15.0;
    gMC->Gsvolu("DCBA", "BOX", idtmed[1809], pbox, 3);
//
// Coil: straight sections, horizontal
    pbox[0] =   6.;
    pbox[1] =  11.;
    pbox[2] = 135.;
    gMC->Gsvolu("DCH1", "BOX", idtmed[1816], pbox, 3);
//
// Coil: straight sections, horizontal
    pbox[0] =   6.;
    pbox[1] =  11.;
    pbox[2] = 135.;
    gMC->Gsvolu("DCH2", "BOX", idtmed[1816], pbox, 3);

//
// Mother volume containing upper coil
    pbox[0] =    8.0;
    pbox[1] =   17.5;
    pbox[2] =  135.0;
    gMC->Gsvolu("DCMU", "BOX", idtmed[1809], pbox, 3);

//
// Coil: straight sections, vertical
    pbox[0] =  6.0;
    pbox[1] =  9.5;
    pbox[2] = 11.0;
    
    gMC->Gsvolu("DCCV", "BOX", idtmed[1816], pbox, 3);
//
// Coil: circular section 

    Float_t ptubs[5];
    ptubs[0] =  0.;
    ptubs[1] = 35.;
    ptubs[2] =  8.;
    ptubs[3] =  0.;
    ptubs[4] = 90.;
//    gMC->Gsvolu("DCC1", "TUBS", idtmed[1809], ptubs, 5);
    ptubs[0] = 13.;
    ptubs[1] = 35.;
    ptubs[2] =  6.;
    ptubs[3] =  0.;
    ptubs[4] = 90.;
    gMC->Gsvolu("DCC1", "TUBS", idtmed[1816], ptubs, 5);
//
// Clamps
    Float_t ppgon[10];
    ppgon[0] =  0.;
    ppgon[1] = 90.;
    ppgon[2] =  1.;
    ppgon[3] =  2.;
    ppgon[4] = -1.;
    ppgon[5] =  0.;
    ppgon[6] = 24.75;
    ppgon[7] =  1.;
    ppgon[8] =  0.;
    ppgon[9] = 24.75;
    gMC->Gsvolu("DCLA", "PGON", idtmed[1809], ppgon, 10);
//
// Assemble all
//
    AliMatrix(idrotm[1811], -90., 0., 90., 90.,   0., 0.);
    AliMatrix(idrotm[1812],   0., 0., 90., 90.,  90., 0.);  
    AliMatrix(idrotm[1813], 180., 0., 90., 90.,  90., 0.);
    AliMatrix(idrotm[1814],   0., 180., 90., 270.,  90., 0.);
    AliMatrix(idrotm[1815], 180., 180., 90., 270.,  90., 0.);  
	
    gMC->Gspos("DCH1", 1, "DCML", 23.25, -13., -17.5, 0, "ONLY");
    gMC->Gspos("DCCV", 1, "DCM0",  12., 19., -159., 0, "ONLY");
    gMC->Gspos("DCCV", 2, "DCM0", -12., 19., -159., 0, "ONLY");
    gMC->Gspos("DCCV", 3, "DCML", 23.25, 20.5, 141.5, 0, "ONLY");

    gMC->Gspos("DCML", 1, "DCM0", -33.25, -2.5, 17.5, 0, "ONLY");
    gMC->Gspos("DCML", 2, "DCM0",  33.25, -2.5, 17.5, idrotm[1811], "ONLY");


    gMC->Gspos("DCH2", 1, "DCMU", 2., 6.5, 0., 0, "ONLY");
    gMC->Gspos("DCMU", 1, "DCM0", -12., 45., 0., 0, "ONLY");
    gMC->Gspos("DCMU", 2, "DCM0",  12., 45., 0., idrotm[1811], "ONLY");

//    gMC->Gspos("DCC2", 1, "DCC1", 0., 0., 0., 0, "ONLY");
    
    gMC->Gspos("DCC1", 1, "DCM0", -12., 27.5,  135., idrotm[1812], "ONLY");
    gMC->Gspos("DCC1", 2, "DCM0",  12., 27.5,  135., idrotm[1812], "ONLY");
    gMC->Gspos("DCC1", 3, "DCM0", -12., 27.5, -135., idrotm[1813], "ONLY");
    gMC->Gspos("DCC1", 4, "DCM0",  12., 27.5, -135., idrotm[1813], "ONLY");

    gMC->Gspos("DCC1", 5, "DCM0",  12., 27.5-32.+13., -135., idrotm[1815], "ONLY");
    gMC->Gspos("DCC1", 6, "DCM0", -12., 27.5-32.+13., -135., idrotm[1815], "ONLY");
    
    gMC->Gspos("DCC1", 7, "DCML", 23.25, -13+13.+11., 117.5, idrotm[1814], "ONLY");

    gMC->Gspos("DCLA", 1, "DCM0",  20., 27.5, -134., 0, "ONLY");
    gMC->Gspos("DCLA", 2, "DCM0",  20., 27.5,  -44., 0, "ONLY");
    gMC->Gspos("DCLA", 3, "DCM0",  20., 27.5,   46., 0, "ONLY");
    gMC->Gspos("DCLA", 4, "DCM0",  20., 27.5,  134., 0, "ONLY");

    gMC->Gspos("DCLA", 5, "DCM0",  -20., 27.5, -134., idrotm[1811], "ONLY");
    gMC->Gspos("DCLA", 6, "DCM0",  -20., 27.5,  -44., idrotm[1811], "ONLY");
    gMC->Gspos("DCLA", 7, "DCM0",  -20., 27.5,   46., idrotm[1811], "ONLY");
    gMC->Gspos("DCLA", 8, "DCM0",  -20., 27.5,  134., idrotm[1811], "ONLY");


    gMC->Gspos("DCBA", 1, "DCM0",  0., -47.5 , 17.5, 0, "ONLY");
    AliMatrix(idrotm[1816], 270., 0., 90., 90.,  180., 0.);  
    gMC->Gspos("DCM0", 1, "ALIC",  0., -6.75,  1080., idrotm[1816], "ONLY");


}

//_____________________________________________________________________________
void AliDIPOv2::DrawModule() const
{
  //
  // Draw a shaded view of the muon absorber
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("DDIP","seen",0);
  gMC->Gsatt("DC1 ","seen",1);
  gMC->Gsatt("DC2 ","seen",1);
  gMC->Gsatt("DC3 ","seen",1);
  gMC->Gsatt("DC4 ","seen",1);
  gMC->Gsatt("DC11","seen",1);
  gMC->Gsatt("DC21","seen",1);
  gMC->Gsatt("DC12","seen",1);
  gMC->Gsatt("DC22","seen",1);
  gMC->Gsatt("DL1 ","seen",1);
  gMC->Gsatt("DL2 ","seen",1);
  gMC->Gsatt("DY1 ","seen",1);
  gMC->Gsatt("DY2 ","seen",1);
  gMC->Gsatt("DYL ","seen",1);
  gMC->Gsatt("DY3 ","seen",1);
 // gMC->Gsatt("DY4 ","seen",1);
 // gMC->Gsatt("DY5 ","seen",1);
 // gMC->Gsatt("DY6 ","seen",1);
//  gMC->Gsatt("DY7 ","seen",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox(".");
  gMC->DefaultRange();
  gMC->Gdraw("alic", 30, 30, 0, 17, 13.5, .019, .019);
  gMC->Gdhead(1111, "Magnetic Dipole Version 2");
  gMC->Gdman(16, 4, "MAN");
}




//_____________________________________________________________________________
void AliDIPOv2::CreateMaterials()
{
  //
  // Create Materials for Magnetic Dipole version 2
  //
  
  Int_t isxfld1   = gAlice->Field()->Integ();
  Int_t isxfld2   = gAlice->Field()->PrecInteg();
  Float_t sxmgmx  = gAlice->Field()->Max();
  
  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };

  Float_t acoil[3]  = { 26.98,1.01,16. };
  Float_t zcoil[3]  = { 13.,1.,8. };
  Float_t wcoil[3]  = { .66,.226,.114 };

  Float_t aresi[3]  = { 1.01,12.011,16.};
  Float_t zresi[3]  = { 1.,6.,8. };
  Float_t wresi[3]  = { .0644,.7655,.1701 };

  Float_t aG10[5] = { 1.01,12.011,16.,28.085 ,79.904 };
  Float_t zG10[5] = { 1.,6.,8.,14.,35. };
  Float_t wG10[5] = { .02089,.22338,.28493,.41342,.05738 };

    // AIR

  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  Float_t dAir1 = 1.20479E-10;

  
  Float_t epsil, stmin, deemax, tmaxfd, stemax;
  
  // --- Define the various materials for GEANT --- 
  //     Aluminum 
  AliMaterial( 9, "ALUMINIUM0$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(29, "ALUMINIUM1$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(49, "ALUMINIUM2$", 26.98, 13., 2.7, 8.9, 37.2);
  
  //     Iron 
  AliMaterial(10, "IRON0$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(30, "IRON1$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(50, "IRON2$     ", 55.85, 26., 7.87, 1.76, 17.1);
  //     Copper
  AliMaterial(17, "COPPER0$   ", 63.55, 29., 8.96, 1.43, 15.1);
  AliMaterial(37, "COPPER1$   ", 63.55, 29., 8.96, 1.43, 15.1);
  AliMaterial(57, "COPPER2$   ", 63.55, 29., 8.96, 1.43, 15.1);
  //     Air 
  AliMixture(15, "AIR0$      ", aAir, zAir, dAir, 4, wAir);
  AliMixture(35, "AIR1$      ", aAir, zAir, dAir, 4, wAir);
  AliMixture(55, "AIR2$      ", aAir, zAir, dAir, 4, wAir);
  AliMixture(75, "AIR_MUON   ", aAir, zAir, dAir, 4, wAir);
  //     Vacuum 
  AliMixture(16, "VACUUM0$ ", aAir, zAir, dAir1, 4, wAir);
  AliMixture(36, "VACUUM1$ ", aAir, zAir, dAir1, 4, wAir);
  AliMixture(56, "VACUUM2$ ", aAir, zAir, dAir1, 4, wAir);
  
  //     stainless Steel 
  AliMixture(19, "STAINLESS STEEL0$", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(39, "STAINLESS STEEL1$", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(59, "STAINLESS STEEL2$", asteel, zsteel, 7.88, 4, wsteel);
  
  //     Coil 
  AliMixture(14, "Al0$", acoil, zcoil, 2.122, 3, wcoil);
  AliMixture(34, "Al1$", acoil, zcoil, 2.122, 3, wcoil);
  AliMixture(54, "Al2$", acoil, zcoil, 2.122, 3, wcoil);

  //RESIN
  AliMixture(13, "RESIN0$", aresi, zresi, 1.05, 3, wresi);
  AliMixture(33, "RESIN1$", aresi, zresi, 1.05, 3, wresi);
  AliMixture(53, "RESIN2$", aresi, zresi, 1.05, 3, wresi);  

  //G10
  AliMixture(11, "G100$", aG10, zG10, 1.7, 5, wG10);
  AliMixture(31, "G101$", aG10, zG10, 1.7, 5, wG10);
  AliMixture(51, "G102$", aG10, zG10, 1.7, 5, wG10); 
 
  // **************** 
  //     Defines tracking media parameters. 
  //     Les valeurs sont commentees pour laisser le defaut 
  //     a GEANT (version 3-21, page CONS200), f.m. 
  epsil  = .001;  // Tracking precision, 
  stemax = -1.;   // Maximum displacement for multiple scat 
  tmaxfd = -20.;  // Maximum angle due to field deflection 
  deemax = -.3;   // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  // *************** 
  
  //    Aluminum 
  AliMedium( 9, "ALU_C0          ",  9, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(29, "ALU_C1          ", 29, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(49, "ALU_C2          ", 49, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Iron 
  AliMedium(10, "FE_C0           ", 10, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(30, "FE_C1           ", 30, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(50, "FE_C2           ", 50, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Air 
  AliMedium(15, "AIR_C0          ", 15, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(35, "AIR_C1          ", 35, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(55, "AIR_C2          ", 55, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(75, "AIR_MUON        ", 75, 0, isxfld2, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Vacuum 
  AliMedium(16, "VA_C0           ", 16, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(36, "VA_C1           ", 36, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(56, "VA_C2           ", 56, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Steel 
  AliMedium(19, "ST_C0           ", 19, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(39, "ST_C1           ", 39, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(59, "ST_C3           ", 59, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Coil 
  AliMedium(14, "Coil_C1         ", 14, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(34, "Coil_C2         ", 34, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(54, "Coil_C3         ", 54, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

  //    Resin 
  AliMedium(13, "RESIN_C0        ", 13, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(33, "RESIN_C1        ", 33, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(53, "RESIN_C2        ", 53, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

  //    G10 
  AliMedium(11, "G10_C0          ", 11, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(31, "G10_C1          ", 31, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(51, "G10_C2          ", 51, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Copper
  AliMedium(17, "Cu_C0           ", 17, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(37, "Cu_C1           ", 37, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(57, "Cu_C2           ", 57, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

}






