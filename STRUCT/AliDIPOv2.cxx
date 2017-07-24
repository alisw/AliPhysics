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

#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>
#include <TGeoVolume.h>
#include <TVirtualMC.h>
#include <TArrayI.h>

#include "AliConst.h"
#include "AliDIPOv2.h"
#include "AliMagF.h"
#include "AliRun.h"

#include "TGeoXtru.h"
#include "TGeoTube.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
 
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
  TVirtualMC::GetMC()->Gsvolu("DDIP", "PCON", idtmed[1874], tpar, 15);
  //
  //       Coils 
  // air - m.f. 
  cpar[0] = 207.;
  cpar[1] = 274.;
  cpar[2] = 37.65;
  cpar[3] = 119.;
  cpar[4] = 241.; 
  //   coil - high cuts
  TVirtualMC::GetMC()->Gsvolu("DC1 ", "TUBS", idtmed[kCoil+40], cpar, 5);
  cpar[3] = -61.;
  cpar[4] = 61.;
  TVirtualMC::GetMC()->Gsvolu("DC2 ", "TUBS", idtmed[kCoil+40], cpar, 5);

  //  coil - low cuts cuts
  cpar[0] = 207.;
//  cpar[1] = cpar[0] + 10.;
  cpar[1] = 217;
  cpar[3] = 119.;
  cpar[4] = 241.;

  TVirtualMC::GetMC()->Gsvolu("DC3 ", "TUBS", idtmed[kCoil], cpar, 5);
  cpar[0] = 207.; 
  cpar[1] = 217; 
  cpar[3] = -61.;
  cpar[4] = 61.;
  TVirtualMC::GetMC()->Gsvolu("DC4 ", "TUBS", idtmed[kCoil], cpar, 5);

  TVirtualMC::GetMC()->Gspos("DC3 ", 1, "DC1 ", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("DC4 ", 1, "DC2 ", 0., 0., 0., 0, "ONLY");

//  dz =  37.65 - 243.55
  dz = -205.9-2.45;
  dx = 5.;
  TVirtualMC::GetMC()->Gspos("DC1 ", 1, "DDIP",  dx, 0.,  dz+kZDipole, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("DC1 ", 2, "DDIP",  dx, 0., -dz+kZDipole, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("DC2 ", 1, "DDIP", -dx, 0.,  dz+kZDipole, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("DC2 ", 2, "DDIP", -dx, 0., -dz+kZDipole, 0, "ONLY");
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
  TVirtualMC::GetMC()->Gsvolu("DC11", "TUBS", idtmed[kCoil+40], cpar, 5);

  dx = TMath::Sin(30.5*kDegrad) * -(207.+33.5)+5./TMath::Sin(30.5*kDegrad); 
  dy = TMath::Cos(30.5*kDegrad) * -(207.+33.5);  
  dz = cpar[1] - 243.55-2.45;
  TVirtualMC::GetMC()->Gspos("DC11", 1, "DDIP",  dx, dy,  dz+kZDipole, idrotm[1800], "ONLY");
  TVirtualMC::GetMC()->Gspos("DC11", 2, "DDIP",  dx, dy, -dz+kZDipole, idrotm[1802], "ONLY");
  TVirtualMC::GetMC()->Gspos("DC11", 3, "DDIP", -dx, dy,  dz+kZDipole, idrotm[1801], "ONLY");
  TVirtualMC::GetMC()->Gspos("DC11", 4, "DDIP", -dx, dy, -dz+kZDipole, idrotm[1803], "ONLY");



//* ... higher cuts
  cpar[0] = 25.;
  cpar[1] = 100.3; //25+75.3
  cpar[2] = 33.5;
  cpar[3] = 0.;
  cpar[4] = 90.;
//*  coil high cuts
  TVirtualMC::GetMC()->Gsvolu("DC12", "TUBS", idtmed[kCoil+40], cpar, 5);

  dx = TMath::Sin(30.5*kDegrad) * -(207.+33.5)+5./TMath::Sin(30.5*kDegrad); 
  dy = TMath::Cos(30.5*kDegrad) *(207.+33.5);  
  dz = cpar[1] - 243.55-2.45;
  TVirtualMC::GetMC()->Gspos("DC12", 1, "DDIP",  dx, dy,  dz+kZDipole, idrotm[1801], "ONLY");
  TVirtualMC::GetMC()->Gspos("DC12", 2, "DDIP",  dx, dy, -dz+kZDipole, idrotm[1803], "ONLY");
  TVirtualMC::GetMC()->Gspos("DC12", 3, "DDIP", -dx, dy,  dz+kZDipole, idrotm[1800], "ONLY");
  TVirtualMC::GetMC()->Gspos("DC12", 4, "DDIP", -dx, dy, -dz+kZDipole, idrotm[1802], "ONLY");

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
  TVirtualMC::GetMC()->Gsvolu("DL1 ", "BOX ", idtmed[kCoil+40], tpar, 3);

// coil - low cuts

  tpar[0] = 5.;
  dx = 37.65  - 5.;  
  TVirtualMC::GetMC()->Gsvolu("DL2 ", "BOX ", idtmed[kCoil], tpar, 3);
  TVirtualMC::GetMC()->Gspos("DL2 ", 1, "DL1 ", dx, 0., 0., 0, "ONLY");

  dx =-53.62;
  dy =-241.26819;
  dz = 0.; 
  TVirtualMC::GetMC()->Gspos("DL1 ", 1, "DDIP", dx,  dy, dz+kZDipole, idrotm[1804], "ONLY");
  TVirtualMC::GetMC()->Gspos("DL1 ", 2, "DDIP", dx, -dy, dz+kZDipole, idrotm[1805], "ONLY");
  TVirtualMC::GetMC()->Gspos("DL1 ", 3, "DDIP",-dx,  dy, dz+kZDipole, idrotm[1806], "ONLY");
  TVirtualMC::GetMC()->Gspos("DL1 ", 4, "DDIP",-dx, -dy, dz+kZDipole, idrotm[1807], "ONLY");

  // Contactor

 //  high cuts

  //Steel outer face planes

  cpar[0] = 207.-18.6;
  cpar[1] = 274.+18.6;
  cpar[2] = 1.;
  cpar[3] = -50.;
  cpar[4] = 50.; 
 
  TVirtualMC::GetMC()->Gsvolu("DCO1", "TUBS", idtmed[1818], cpar, 5);
  dx = -5.;
  dz = 168.25-1.5-1.;
  TVirtualMC::GetMC()->Gspos("DCO1", 1, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");
  dz = 243.55+4.5+1.5+1.;
  TVirtualMC::GetMC()->Gspos("DCO1", 2, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");
  
  // 9.06.2000

  //  cpar[0] = 207.-18.6;
  //  cpar[1] = 274.+18.6;
  // cpar[2] = 1.;
  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 
 
  TVirtualMC::GetMC()->Gsvolu("DCO2", "TUBS", idtmed[1818], cpar, 5);
  dx = +5.;
  dz = 168.25-1.5-1.;
  TVirtualMC::GetMC()->Gspos("DCO2", 1, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");
  dz = 243.55+4.5+1.5+1.;
  TVirtualMC::GetMC()->Gspos("DCO2", 2, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");
 


  // Resin face planes

  cpar[0] = 207.;
  cpar[1] = 274.;
  cpar[2] = .75;
  cpar[3] = -50.;
  cpar[4] = 50.; 
 
  TVirtualMC::GetMC()->Gsvolu("DCO3", "TUBS", idtmed[1812], cpar, 5);
  dx = -5;
  dz = 168.25-0.75;
  TVirtualMC::GetMC()->Gspos("DCO3", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");
  dz = 243.55+4.5+0.75;
  TVirtualMC::GetMC()->Gspos("DCO3", 2, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000

  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 
  TVirtualMC::GetMC()->Gsvolu("DCO4", "TUBS", idtmed[1812], cpar, 5);
  dx = +5;
  dz = 168.25-0.75;
  TVirtualMC::GetMC()->Gspos("DCO4", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");
  dz = 243.55+4.5+0.75 ;
  TVirtualMC::GetMC()->Gspos("DCO4", 2, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

 
  // G10 face plane

  cpar[0] = 207.;
  cpar[1] = 274.;
  cpar[2] = 2.25;
  cpar[3] = -50.;
  cpar[4] = 50.; 
 
  TVirtualMC::GetMC()->Gsvolu("DCO5", "TUBS", idtmed[1810], cpar, 5);

  dx = -5;
  dz = 243.55+2.25;
  TVirtualMC::GetMC()->Gspos("DCO5", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000

  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 

  TVirtualMC::GetMC()->Gsvolu("DCO6", "TUBS", idtmed[1810], cpar, 5);

  dx = +5;
  dz = 243.55+2.25;
  TVirtualMC::GetMC()->Gspos("DCO6", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  //Steel supported planes

  cpar[0] = 274.+1.5+2.;
  cpar[1] = 274.+18.6;
  cpar[2] = 1.;
  cpar[3] = -50.;
  cpar[4] = 50.;  
 
  TVirtualMC::GetMC()->Gsvolu("DCO7", "TUBS", idtmed[1818], cpar, 5);

  dx = -5;
  dz = 168.25+1.;
  TVirtualMC::GetMC()->Gspos("DCO7", 1, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000
  cpar[0] = 274.+1.5+2.;
  cpar[1] = 274.+18.6;
  cpar[2] = 1.;
  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 

 
  TVirtualMC::GetMC()->Gsvolu("DCO8", "TUBS", idtmed[1818], cpar, 5);

  dx = +5;
  dz = 168.25+1.;
  TVirtualMC::GetMC()->Gspos("DCO8", 1, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");

  //

  cpar[0] = 207.- 18.6;
  cpar[1] = 207.- 2.- 1.5;
  cpar[2] = 1.;
  cpar[3] = -50.;
  cpar[4] = 50.; 

  TVirtualMC::GetMC()->Gsvolu("DCO9", "TUBS", idtmed[1818], cpar, 5);

  dx = -5;
  dz = 168.25+1.;
  TVirtualMC::GetMC()->Gspos("DCO9", 1, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000

  cpar[0] = 207.- 18.6;
  cpar[1] = 207.- 2.- 1.5;
  cpar[2] = 1.;
  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 

  TVirtualMC::GetMC()->Gsvolu("DCOA", "TUBS", idtmed[1818], cpar, 5);

  dx = +5;
  dz = 168.25+1.;
  TVirtualMC::GetMC()->Gspos("DCOA", 1, "DDIP", dx, 0, dz+kZDipole, 0, "ONLY");


  // Sides steel planes

  cpar[0] = 207. - 1.5 -2.;
  cpar[1] = 207. - 1.5;
  cpar[2] = ((243.55+4.5+1.5)-168.25)/2;
  cpar[3] = -50.;
  cpar[4] = 50.; 
 
  TVirtualMC::GetMC()->Gsvolu("DCOB", "TUBS", idtmed[1818], cpar, 5);

  cpar[0] = 274. + 1.5;
  cpar[1] = 274. + 1.5 +2.;

  TVirtualMC::GetMC()->Gsvolu("DCOC", "TUBS", idtmed[1818], cpar, 5);

  dx=-5.;
  dz = ((243.55+4.5+1.5)+168.25)/2;
  TVirtualMC::GetMC()->Gspos("DCOB", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("DCOC", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000

  cpar[0] = 207. - 1.5 -2.;
  cpar[1] = 207. - 1.5;
  cpar[2] = ((243.55+4.5+1.5)-168.25)/2;
  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 

  TVirtualMC::GetMC()->Gsvolu("DCOD", "TUBS", idtmed[1818], cpar, 5);

  cpar[0] = 274. + 1.5;
  cpar[1] = 274. + 1.5 +2.;

  TVirtualMC::GetMC()->Gsvolu("DCOE", "TUBS", idtmed[1818], cpar, 5);

  dx=+5.;
  dz = ((243.55+4.5+1.5)+168.25)/2;
  TVirtualMC::GetMC()->Gspos("DCOD", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("DCOE", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");


  // Top and bottom resin  planes

  cpar[0] = 207. - 1.5;
  cpar[1] = 207.;
  cpar[2] = ((243.55+4.5+1.5)-168.25)/2;
  cpar[3] = -50.;
  cpar[4] = 50.; 
 
  TVirtualMC::GetMC()->Gsvolu("DCOF", "TUBS", idtmed[1812], cpar, 5);

  cpar[0] = 274.;
  cpar[1] = 274. + 1.5;

  TVirtualMC::GetMC()->Gsvolu("DCOG", "TUBS", idtmed[1812], cpar, 5);


  dx=-5.;
  dz = ((243.55+4.5+1.5)+168.25)/2;
  TVirtualMC::GetMC()->Gspos("DCOF", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("DCOG", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000
  cpar[0] = 207. - 1.5;
  cpar[1] = 207.;
  cpar[2] = ((243.55+4.5+1.5)-168.25)/2;

  cpar[3] = 180.-50.;
  cpar[4] = 180.+50.; 

  TVirtualMC::GetMC()->Gsvolu("DCOH", "TUBS", idtmed[1812], cpar, 5);

  cpar[0] = 274.;
  cpar[1] = 274. + 1.5;

  TVirtualMC::GetMC()->Gsvolu("DCOI", "TUBS", idtmed[1812], cpar, 5);


  dx=+5.;
  dz = ((243.55+4.5+1.5)+168.25)/2;
  TVirtualMC::GetMC()->Gspos("DCOH", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("DCOI", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");


  // Aluminum cabels

  cpar[0] = 274. + 1.5  +2.;
  cpar[1] = 274. + 1.5  +2. + 80.;
  cpar[2] = 5.05/2;
  cpar[3] = -24.;
  cpar[4] = 24.; 
 
  TVirtualMC::GetMC()->Gsvolu("DCOJ", "TUBS", idtmed[kCable], cpar, 5);

  //  dx = 274. + 1.5  +2. +40.;
  //  dx = 5. + 1.5 +2. +40.;
  //  dx = 5. + 1.5 +2.;
  dx=-5.;
  dz = 168.25 + 5.05 + 5.05/2;
  TVirtualMC::GetMC()->Gspos("DCOJ", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  dz = 243.55 - 5.05/2;
  TVirtualMC::GetMC()->Gspos("DCOJ", 2, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  // 9.06.2000

  cpar[3] = 180.-24.;
  cpar[4] = 180.+24.; 

   TVirtualMC::GetMC()->Gsvolu("DCOK", "TUBS", idtmed[kCable], cpar, 5);

  //  dx = 274. + 1.5  +2. +40.;
  //  dx = 5. + 1.5 +2. +40.;
  //  dx = 5. + 1.5 +2.;
  dx=+5.;
  dz = 168.25 + 5.05 + 5.05/2;
  TVirtualMC::GetMC()->Gspos("DCOK", 1, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

  dz = 243.55 - 5.05/2;
  TVirtualMC::GetMC()->Gspos("DCOK", 2, "DDIP", dx,  0, dz+kZDipole, 0, "ONLY");

 
  //   YOKE 

// Top and bottom blocks
  ypar[0] = 298.1; 
  ypar[1] = 69.5;
  ypar[2] = 155.75;

// iron- high cuts
  TVirtualMC::GetMC()->Gsvolu("DY1 ", "BOX ", idtmed[1858], ypar, 3);
  ypar[0] = 144.+10.; 
  ypar[1] = 193.3+10.;
  ypar[2] = 5.;
  ypar[3] = 155.75;
  dy = -69.5 + 5.;
// iron- low cuts
  TVirtualMC::GetMC()->Gsvolu("DY11", "TRD1", idtmed[1818], ypar, 4);
  TVirtualMC::GetMC()->Gspos("DY11", 1, "DY1 ", 0.,  dy, 0., 0, "ONLY");

  dy = 365.5;
  dz = 4.95;
  TVirtualMC::GetMC()->Gspos("DY1 ", 1, "DDIP", 0.,  dy, -dz+kZDipole, 0, "ONLY");

  the1 = 270.;
  phi1 = 0.;
  the2 = 270.;
  phi2 = 90.;
  the3 = 0.;
  phi3 = 0.;
  AliMatrix(idrotm[1808], the1, phi1, the2, phi2, the3, phi3);
  TVirtualMC::GetMC()->Gspos("DY1 ", 2, "DDIP", 0., -dy, -dz+kZDipole, idrotm[1808] , "ONLY");

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

  TVirtualMC::GetMC()->Gsvolu("DY2 ", "TRAP", idtmed[1858], ypar,11);

  ypar[4] = 47.9 -5.;
  ypar[5] = 72.55 -5.;

  ypar[8] = 47.9 -5.;
  ypar[9] = 72.55 -5.;


// iron - low cuts

  TVirtualMC::GetMC()->Gsvolu("DY22", "TRAP", idtmed[1818], ypar,11);

  dy = 0.;
  dx = -5.;

  TVirtualMC::GetMC()->Gspos("DY22", 1, "DY2 ", dx,  dy, 0., 0, "ONLY");

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
  
  TVirtualMC::GetMC()->Gspos("DY2 ", 1, "DDIP",  dx, 0.0,  dz+kZDipole, idrotm[1809], "ONLY");
  TVirtualMC::GetMC()->Gspos("DY2 ", 2, "DDIP", -dx, 0.0,  dz+kZDipole, idrotm[1810], "ONLY");
  
  AliMatrix(idrotm[1811], 270., 0., 90., 90., 180., 0.);
  TVirtualMC::GetMC()->Gspos("DDIP", 1, "ALIC", 0., 0., 0., idrotm[1811], "ONLY");
  gGeoManager->SetVolumeAttribute("DDIP", "SEEN", 0);
}


void AliDIPOv2::CreateCompensatorDipole()
{
    //
    //  Geometry of the compensator Dipole MBWMD (was MCB @ SPS)
    //  LAB I/EA Note 74.10
    //  6/5/1974
    //
    const Float_t kHCoil       =  22.;  // Coil Height
    const Float_t kWCoil       =  12.;  // Coil Width
    const Float_t kLCoilH      = 250.;  // Hor. Coil Length
    const Float_t kRCoilC      =  31.;  // Circ Coil Radius
    const Float_t kWBase       = 125.;  // Base Width
    const Float_t kHBase       =  30.;  // Base Height
    
    const Float_t kWUYoke      =  16.;
    const Float_t kHUYoke      =  31.;

    const Float_t kWLYoke      =  50.0;
    const Float_t kHLYoke      =  61.0;
    const Float_t kLLYoke      =  kLCoilH + kRCoilC;

    const Float_t kWApperture  =  12.;
    const Float_t kDCoil       =  kHUYoke + kHLYoke - 6. - 2. * kRCoilC;
    
    const Float_t kH           =  kHBase +  kHUYoke +  kHLYoke;
    
    Int_t *idtmed = fIdtmed->GetArray()-1799;
    Int_t idrotm[1899];
//
    Float_t pbox[3];
//  Mother volumes
/*    
    goto ECV_CODE;
    TGeoVolumeAssembly* asDCM0 = new TGeoVolumeAssembly("DCM0");
    asDCM0->SetName("DCM0");
//
//  Mother volume containing lower coil
    pbox[0] = kWLYoke / 2.;
    pbox[1] = kHLYoke / 2.;
    pbox[2] = kLLYoke / 2.;
    
    TVirtualMC::GetMC()->Gsvolu("DCML", "BOX", idtmed[1809 + 40], pbox, 3);
//
// Base
    pbox[0] = kWBase / 2.;
    pbox[1] = kHBase / 2.;
    TVirtualMC::GetMC()->Gsvolu("DCBA", "BOX", idtmed[1809 + 40], pbox, 3);
//
// Coil: straight sections, horizontal
    pbox[0] = kWCoil  / 2.;
    pbox[1] = kHCoil  / 2.;
    pbox[2] = kLCoilH / 2.;
    TVirtualMC::GetMC()->Gsvolu("DCH1", "BOX", idtmed[1816 + 40], pbox, 3);
//
// Coil: straight sections, horizontal
    pbox[0] = kWCoil  / 2.;
    pbox[1] = kHCoil  / 2.;
    pbox[2] = kLCoilH / 2.;
    TVirtualMC::GetMC()->Gsvolu("DCH2", "BOX", idtmed[1816 + 40], pbox, 3);

//
// Mother volume containing upper coil
    pbox[0] =  kWUYoke / 2.;
    pbox[1] =  kHUYoke / 2.;
    pbox[2] =  kLCoilH / 2.;
    TVirtualMC::GetMC()->Gsvolu("DCMU", "BOX", idtmed[1809 + 40], pbox, 3);

//
// Coil: straight sections, vertical
    pbox[0] = kWCoil / 2.;
    pbox[1] = kDCoil / 2.;
    pbox[2] = kHCoil / 2.;
    
    TVirtualMC::GetMC()->Gsvolu("DCCV", "BOX", idtmed[1816 + 40], pbox, 3);
//
// Coil: circular section 

    Float_t ptubs[5];
    ptubs[0] = kRCoilC - kHCoil;
    ptubs[1] = kRCoilC;
    ptubs[2] = kWCoil / 2.;
    ptubs[3] =  0.;
    ptubs[4] = 90.;
    TVirtualMC::GetMC()->Gsvolu("DCC1", "TUBS", idtmed[1816 + 40], ptubs, 5);
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
    TVirtualMC::GetMC()->Gsvolu("DCLA", "PGON", idtmed[1809 + 40], ppgon, 10);
//
// Assemble all
//
    AliMatrix(idrotm[1811], -90., 0., 90., 90.,   0., 0.);
    AliMatrix(idrotm[1812],   0., 0., 90., 90.,  90., 0.);  
    AliMatrix(idrotm[1813], 180., 0., 90., 90.,  90., 0.);
    AliMatrix(idrotm[1814],   0., 180., 90., 270.,  90., 0.);
    AliMatrix(idrotm[1815], 180., 180., 90., 270.,  90., 0.);  
    
    Float_t dx, dy, dz;
    Float_t dy0 = 0.;

    dy0 = -kH / 2. + kHBase/2.;
    TVirtualMC::GetMC()->Gspos("DCBA", 1, "DCM0",  0., dy0, 15.0, 0, "ONLY");
    
    // Lower coil
    dx = ( kWLYoke - kWCoil) / 2.;
    dy = (-kHLYoke + kHCoil) / 2. + 6.;
    TVirtualMC::GetMC()->Gspos("DCH1", 1, "DCML",  dx, dy,  -kRCoilC / 2., 0, "ONLY");
    // Lower mother volume
    dx = (kWLYoke + kWApperture) / 2.;
    dy0 += (kHBase +  kHLYoke) / 2.;
    TVirtualMC::GetMC()->Gspos("DCML", 1, "DCM0", -dx, dy0, kRCoilC / 2., 0, "ONLY");
    TVirtualMC::GetMC()->Gspos("DCML", 2, "DCM0", +dx, dy0, kRCoilC / 2., idrotm[1811], "ONLY");
    
    dx = (kWUYoke - kWCoil) / 2.;
    dy = (kHUYoke - kHCoil) / 2;
    // Upper coil
    TVirtualMC::GetMC()->Gspos("DCH2", 1, "DCMU",   dx,  dy, 0., 0, "ONLY");
    // Upper mother volume
    dx = (kWUYoke + kWApperture) / 2.;
    dy0 +=  (kHLYoke + kHUYoke) / 2.;
    TVirtualMC::GetMC()->Gspos("DCMU", 1, "DCM0", -dx, dy0, 0., 0, "ONLY");
    TVirtualMC::GetMC()->Gspos("DCMU", 2, "DCM0", +dx, dy0, 0., idrotm[1811], "ONLY");

    // Vertical coils
    dx =  (kWCoil +  kWApperture) / 2.;
    dy =  kH / 2. - kDCoil / 2. - kRCoilC;
    dz =  (kLCoilH - kHCoil) / 2. + kRCoilC;
    TVirtualMC::GetMC()->Gspos("DCCV", 1, "DCM0",  dx,  dy, -dz, 0, "ONLY");
    TVirtualMC::GetMC()->Gspos("DCCV", 2, "DCM0", -dx,  dy, -dz, 0, "ONLY");

    dx = (kWLYoke - kWCoil) / 2.;
    dy = -kHLYoke / 2. + kDCoil / 2. + 6. + kRCoilC;
    dz =  kLLYoke / 2. - kHCoil / 2.;
    
    TVirtualMC::GetMC()->Gspos("DCCV", 3, "DCML", dx, dy,  dz, 0, "ONLY");



    // Circular coil
    dx =  (kWCoil +  kWApperture) / 2.;
    dy = dy0 + kHUYoke / 2. - kRCoilC;
    dz =  kLCoilH / 2.;
    TVirtualMC::GetMC()->Gspos("DCC1", 1, "DCM0", -dx, dy,  dz, idrotm[1812], "ONLY");
    TVirtualMC::GetMC()->Gspos("DCC1", 2, "DCM0", +dx, dy,  dz, idrotm[1812], "ONLY");
    TVirtualMC::GetMC()->Gspos("DCC1", 3, "DCM0", +dx, dy, -dz, idrotm[1813], "ONLY");
    TVirtualMC::GetMC()->Gspos("DCC1", 4, "DCM0", -dx, dy, -dz, idrotm[1813], "ONLY");
    dy = -kH / 2. + kHBase + 6. + kRCoilC;
    TVirtualMC::GetMC()->Gspos("DCC1", 5, "DCM0", +dx, dy, -dz, idrotm[1815], "ONLY");
    TVirtualMC::GetMC()->Gspos("DCC1", 6, "DCM0", -dx, dy, -dz, idrotm[1815], "ONLY");

    dx = ( kWLYoke - kWCoil) / 2.;
    dy =  -kHLYoke / 2. + 6. + kRCoilC;
    dz =  kLLYoke / 2. - kRCoilC;
    TVirtualMC::GetMC()->Gspos("DCC1", 7, "DCML", dx, dy, dz, idrotm[1814], "ONLY");

//  Clamps
    dx = kWApperture / 2. + kWUYoke;
    dy = -kH / 2. + kHLYoke + kHBase;
    

    TVirtualMC::GetMC()->Gspos("DCLA", 1, "DCM0",  dx, dy, -119., 0, "ONLY");
    TVirtualMC::GetMC()->Gspos("DCLA", 2, "DCM0",  dx, dy,  -44., 0, "ONLY");
    TVirtualMC::GetMC()->Gspos("DCLA", 3, "DCM0",  dx, dy,   46., 0, "ONLY");
    TVirtualMC::GetMC()->Gspos("DCLA", 4, "DCM0",  dx, dy,  119., 0, "ONLY");

    TVirtualMC::GetMC()->Gspos("DCLA", 5, "DCM0",  -dx, dy, -119., idrotm[1811], "ONLY");
    TVirtualMC::GetMC()->Gspos("DCLA", 6, "DCM0",  -dx, dy,  -44., idrotm[1811], "ONLY");
    TVirtualMC::GetMC()->Gspos("DCLA", 7, "DCM0",  -dx, dy,   46., idrotm[1811], "ONLY");
    TVirtualMC::GetMC()->Gspos("DCLA", 8, "DCM0",  -dx, dy,  119., idrotm[1811], "ONLY");


    AliMatrix(idrotm[1816], 270., 0., 90., 90.,  180., 0.);  
    // TVirtualMC::GetMC()->Gspos("DCM0", 1, "ALIC",  0., -12.,  1075., idrotm[1816], "ONLY");

*/
    // ECV:
ECV_CODE: 
    // asDCM0->SetName("DCM0");
    TGeoVolume * ALIC = gGeoManager->GetVolume("ALIC");
    ALIC->AddNode(CreateMagnetYoke(), 1, new TGeoTranslation(0., 0., 1075.));
    //
}

TGeoVolume * AliDIPOv2::CreateMagnetYoke()
{
  TGeoVolumeAssembly * voMagnet = new TGeoVolumeAssembly("DCM0");
  voMagnet->SetName("DCM0");
  TGeoRotation * Ry180 = new TGeoRotation("Ry180", 180., 180.,   0.);  
  TGeoMedium * kMedAlu     = gGeoManager->GetMedium("DIPO_ALU_C0");
  TGeoMedium * kMedCooper  = gGeoManager->GetMedium("DIPO_Cu_C0");
  TGeoMedium * kMedIron    = gGeoManager->GetMedium("DIPO_FE_C0");

  new TGeoBBox("shMagnetYokeOuter"   , 116.4/2.0 , 90.2/2.0 , 250.0/2.0 ); 
  new TGeoBBox("shMagnetYokeInnerUp" ,   8.0/2.0 , 32.2/2.0 , 250.0/1.0  ); 
  new TGeoBBox("shMagnetYokeInnerDw" ,  46.0/2.0 , 23.0/2.0 , 250.0/1.0  ); 
  (new TGeoTranslation("trMagnetYokeOuter"  ,  0.0, -29.1, 0.0)) -> RegisterYourself();
  (new TGeoTranslation("trMagnetYokeInnerUp",  0.0,   0.0, 0.0)) -> RegisterYourself();
  (new TGeoTranslation("trMagnetYokeInnerDw",  0.0, -27.5, 0.0)) -> RegisterYourself();

  TGeoCompositeShape * shMagnetYoke = new TGeoCompositeShape("shMagnet", "shMagnetYokeOuter:trMagnetYokeOuter-(shMagnetYokeInnerUp:trMagnetYokeInnerUp+shMagnetYokeInnerDw:trMagnetYokeInnerDw)");
  TGeoVolume * voMagnetYoke = new TGeoVolume("voMagnetYoke", shMagnetYoke, kMedIron);

  // Make the coils:
  TGeoVolume * voCoilH = gGeoManager -> MakeBox("voCoilH", kMedCooper, 12.64/2.0,  21.46/2.0, 310.5/2.0);
  TGeoVolume * voCoilV = gGeoManager -> MakeBox("voCoilV", kMedCooper, 12.64/2.0,  35.80/2.0,  26.9/2.0);

  // Make the top coil supports:
  // Polygone Coordinates (x,y)
  Double_t x, y;
  const Double_t kDegToRad  = TMath::Pi()/180. ; 
  const Double_t AngleInner = 4.5 * kDegToRad  ; 
  const Double_t AngleOuter = 56.0 * kDegToRad ; 
  const Double_t ArcStart   = 90. - AngleOuter / kDegToRad ; 
  const Double_t ArcEnd     = 90. + AngleInner / kDegToRad ; 
  const Double_t b          = 13.6             ; 
  const Double_t Lx         = 37.2             ; 
  const Double_t Ly         = 25.7             ; 
  const Double_t LxV        = 14.9; 
  const Double_t R          = 9.50;
  const Double_t dz         = 2.00/2.0;
  const Int_t npoints = 8;
  Double_t CenterX; Double_t CenterY; 
  Double_t PointsX[npoints] = {0.};
  Double_t PointsY[npoints] = {0.};
  Int_t ip = 0;
  // Start point: 
  x = 0.0; y = 0.0;
  PointsX[ip]=x; PointsY[ip]=y; ip++;
  // 1st step: 
  x = 0.00;
  y = 1.95;
  PointsX[ip]=x; PointsY[ip]=y; ip++;
  // 2nd step: 
  x+= b;
  y+= b * TMath::Tan(AngleInner);
  PointsX[ip]=x; PointsY[ip]=y; ip++;
  // Center of Arc:
  x+= R * TMath::Sin(AngleInner);
  y-= R * TMath::Cos(AngleInner);
  CenterX=x; CenterY=y; 
  TGeoTubeSeg * shPolygonArc = new TGeoTubeSeg("shPolygonArc", R-2.0, R, dz , ArcStart, ArcEnd);
  (new TGeoTranslation("trPolygonArc", x,y,0.))->RegisterYourself();
  // 3rd Step:
  x+= R * TMath::Sin(AngleOuter);
  y+= R * TMath::Cos(AngleOuter);
  PointsX[ip]=x; PointsY[ip]=y; ip++;
  // 4th Step:
  Double_t a = Lx  - b - R * TMath::Sin(AngleInner) - R * TMath::Sin(AngleOuter);
  x = Lx; 
  y-= a * TMath::Tan(AngleOuter);
  PointsX[ip]=x; PointsY[ip]=y; ip++;
  // 5th Step:
  x =  Lx;
  y = -Ly;
  PointsX[ip]=x; PointsY[ip]=y; ip++;
  // 6th Step:
  x =   LxV;
  y =   -Ly;
  PointsX[ip]=x; PointsY[ip]=y; ip++;
  // 7th Step:
  x =   LxV;
  y =   0.0;
  PointsX[ip]=x; PointsY[ip]=y; ip++;
  // 
  //
  //
  TGeoXtru * shPolygon = new TGeoXtru(2);
  shPolygon->SetNameTitle("shPolygon","shPolygon");
  shPolygon->DefinePolygon(npoints,PointsX,PointsY);
  shPolygon->DefineSection(0, -dz, 0., 0., 1.0); // index, Z position, offset (x,y) and scale for first section
  shPolygon->DefineSection(1, +dz, 0., 0., 1.0); // idem, second section

  TGeoCompositeShape * shCoilSupportV = new TGeoCompositeShape("shCoilSupportV", "shPolygon+shPolygonArc:trPolygonArc");
  TGeoVolume * voCoilSupportV = new TGeoVolume("voCoilSupportV", shCoilSupportV, kMedAlu);
  
  const Double_t MagCoilDx = 12.64/2.;
  const Double_t MagCoilDy = 21.46/2.;
  const Double_t SqOuterDx = MagCoilDx + 2.8;
  const Double_t SqInnerDx = MagCoilDx + 0.6;
  const Double_t SqOuterDy =  29.2/2.;
  const Double_t SqInnerDy =  24.8/2.;
  const Double_t SqOuterDz =  15.5/2.;
  const Double_t SqInnerDz = SqOuterDz * 2.;
  TGeoBBox * shCoilSupportSqOuter = new TGeoBBox("shCoilSupportSqOuter", SqOuterDx, SqOuterDy, SqOuterDz);
  TGeoBBox * shCoilSupportSqInner = new TGeoBBox("shCoilSupportSqInner", SqInnerDx, SqInnerDy, SqInnerDz);
  TGeoCompositeShape * shCoilSupportSq = new TGeoCompositeShape("shCoilSupportSq", "shCoilSupportSqOuter - shCoilSupportSqInner");
  TGeoVolume * voCoilSupportSq = new TGeoVolume("voCoilSupportSq", shCoilSupportSq, kMedAlu);

  const Double_t HSuppDx = (Lx - LxV + 0.6)/2.0;
  const Double_t HSuppDy = 2.2/2.0;
  const Double_t HSuppDz = SqOuterDz; 

  TGeoVolume * voCoilSupportH = gGeoManager -> MakeBox("voCoilSupportH", kMedAlu, HSuppDx, HSuppDy, HSuppDz);
  
  TGeoVolumeAssembly * voCoilSupport = new TGeoVolumeAssembly("voCoilSupport");
  voCoilSupportV  -> SetLineColor(kViolet+9);
  voCoilSupportSq -> SetLineColor(kBlue-5);
  voCoilSupportH  -> SetLineColor(kPink);
  // voCoilSupportH  -> SetTransparency(16);
  voCoilSupport->AddNode(voCoilSupportV  , 1 , new TGeoTranslation(SqOuterDx - LxV     , SqOuterDy                , 0. )); 
  voCoilSupport->AddNode(voCoilSupportSq , 1 , new TGeoTranslation(0.                  , 0.                       , 0. )); 
  voCoilSupport->AddNode(voCoilSupportH  , 1 , new TGeoTranslation(SqOuterDx + HSuppDx , SqOuterDy - Ly - HSuppDy , 0. )); 

  // Make the Top Support for Geodesic reference points:
  TGeoVolume * voSupportHTop = gGeoManager -> MakeBox("voSupportHTop", kMedAlu, 66.0/2.0,   2.0/2.0, 17.0/2.0);
  TGeoVolume * voSupportHBot = gGeoManager -> MakeBox("voSupportHBot", kMedAlu, 14.0/2.0,   2.0/2.0, 17.0/2.0);
  TGeoVolume * voSupportVert = gGeoManager -> MakeBox("voSupportVert", kMedAlu,  3.0/2.0,  25.0/2.0, 17.0/2.0);

  TGeoVolumeAssembly * voSupportGeoRefPoint = new TGeoVolumeAssembly("voSupportGeoRefPoint");
  voSupportHTop -> SetLineColor(kGreen);
  voSupportHBot -> SetLineColor(kGreen);
  voSupportVert -> SetLineColor(kGreen);
  voSupportGeoRefPoint -> AddNode( voSupportHTop , 1 , new TGeoTranslation(  0.0 , 28.0 , 0. )); 
  voSupportGeoRefPoint -> AddNode( voSupportHBot , 1 , new TGeoTranslation(+33.0 ,  1.0 , 0. )); 
  voSupportGeoRefPoint -> AddNode( voSupportHBot , 2 , new TGeoTranslation(-33.0 ,  1.0 , 0. )); 
  voSupportGeoRefPoint -> AddNode( voSupportVert , 1 , new TGeoTranslation(+31.5 , 14.5 , 0. )); 
  voSupportGeoRefPoint -> AddNode( voSupportVert , 2 , new TGeoTranslation(-31.5 , 14.5 , 0. )); 

  // Add some color:
  voMagnetYoke  -> SetLineColor(kAzure-7);
  voCoilH       -> SetLineColor(kOrange-3);
  voCoilV       -> SetLineColor(kOrange-3);
  // Assembling:
  
  voMagnet -> AddNode(voMagnetYoke         , 1 , new TGeoTranslation(0.     , 0.     ,    0.0 )); 
  voMagnet -> AddNode(voCoilH              , 1 , new TGeoTranslation(+16.14 , +29.83 ,    0.0 )); 
  voMagnet -> AddNode(voCoilH              , 2 , new TGeoTranslation(-16.14 , +29.83 ,    0.0 )); 
  voMagnet -> AddNode(voCoilH              , 3 , new TGeoTranslation(+16.14 , -27.43 ,    0.0 )); 
  voMagnet -> AddNode(voCoilH              , 4 , new TGeoTranslation(-16.14 , -27.43 ,    0.0 )); 
  voMagnet -> AddNode(voCoilV              , 1 , new TGeoTranslation(+16.14 , 1.20   , +141.8 )); 
  voMagnet -> AddNode(voCoilV              , 2 , new TGeoTranslation(-16.14 , 1.20   , +141.8 )); 
  voMagnet -> AddNode(voCoilV              , 3 , new TGeoTranslation(+16.14 , 1.20   , -141.8 )); 
  voMagnet -> AddNode(voCoilV              , 4 , new TGeoTranslation(-16.14 , 1.20   , -141.8 )); 
  Double_t zGeoRef = 74.0/2. + SqOuterDz + 9.0 + 17.0/2.0;
  voMagnet -> AddNode(voSupportGeoRefPoint , 1 , new TGeoTranslation( 0.    , 16.0   , +zGeoRef)); 
  voMagnet -> AddNode(voSupportGeoRefPoint , 2 , new TGeoTranslation( 0.    , 16.0   , -zGeoRef)); 
  Double_t zCoilSupp  = 29.83 - MagCoilDy -0.6 + SqInnerDy;
  voMagnet -> AddNode(voCoilSupport        , 1 , new TGeoTranslation(+16.14 , zCoilSupp ,   74.0*0.5 )); 
  voMagnet -> AddNode(voCoilSupport        , 2 , new TGeoTranslation(+16.14 , zCoilSupp ,  -74.0*0.5 )); 
  voMagnet -> AddNode(voCoilSupport        , 3 , new TGeoTranslation(+16.14 , zCoilSupp ,   74.0*1.5 )); 
  voMagnet -> AddNode(voCoilSupport        , 4 , new TGeoTranslation(+16.14 , zCoilSupp ,  -74.0*1.5 )); 
  // 
  voMagnet -> AddNode(voCoilSupport        , 5 , new TGeoCombiTrans (-16.14 , zCoilSupp ,   74.0*0.5, Ry180)); 
  voMagnet -> AddNode(voCoilSupport        , 6 , new TGeoCombiTrans (-16.14 , zCoilSupp ,  -74.0*0.5, Ry180)); 
  voMagnet -> AddNode(voCoilSupport        , 7 , new TGeoCombiTrans (-16.14 , zCoilSupp ,   74.0*1.5, Ry180)); 
  voMagnet -> AddNode(voCoilSupport        , 8 , new TGeoCombiTrans (-16.14 , zCoilSupp ,  -74.0*1.5, Ry180)); 

  return (TGeoVolume*) voMagnet;

/*
  // TGeoVolumeAssembly * voMagnet = new TGeoVolumeAssembly("voMagnet");
  new TGeoBBox("shMagnetYokeOuter"   , 116.4/2.0 , 90.2/2.0 , 250.0/2.0 ); 
  new TGeoBBox("shMagnetYokeInnerUp" ,   8.0/2.0 , 32.2/2.0 , 250.0/1.0  ); 
  new TGeoBBox("shMagnetYokeInnerDw" ,  46.0/2.0 , 23.0/2.0 , 250.0/1.0  ); 
  (new TGeoTranslation("trMagnetYokeOuter"  ,  0.0, -29.1, 0.0)) -> RegisterYourself();
  (new TGeoTranslation("trMagnetYokeInnerUp",  0.0,   0.0, 0.0)) -> RegisterYourself();
  (new TGeoTranslation("trMagnetYokeInnerDw",  0.0, -27.5, 0.0)) -> RegisterYourself();

  TGeoCompositeShape * shMagnetYoke = new TGeoCompositeShape("shMagnet", "shMagnetYokeOuter:trMagnetYokeOuter-(shMagnetYokeInnerUp:trMagnetYokeInnerUp+shMagnetYokeInnerDw:trMagnetYokeInnerDw)");
  TGeoVolume * voMagnetYoke = new TGeoVolume("voMagnetYoke", shMagnetYoke, kMedIron   );

  // Make the coils:
  TGeoVolume * voCoilH = gGeoManager -> MakeBox("voCoilH", kMedCooper, 12.64/2.0,  21.46/2.0, 310.5/2.0);
  TGeoVolume * voCoilV = gGeoManager -> MakeBox("voCoilV", kMedCooper, 12.64/2.0,  35.80/2.0,  26.9/2.0);

  // Add some color:
  voMagnetYoke -> SetLineColor(kAzure-7);
  voCoilH      -> SetLineColor(kOrange-3);
  voCoilV      -> SetLineColor(kOrange-3);
  // Assembling:
  
  voMagnet -> AddNode(voMagnetYoke, 1, new TGeoTranslation(0., 0., 0.));
  voMagnet -> AddNode(voCoilH,      1, new TGeoTranslation(+16.14,+29.83,   0.0));
  voMagnet -> AddNode(voCoilH,      2, new TGeoTranslation(-16.14,+29.83,   0.0));
  voMagnet -> AddNode(voCoilH,      3, new TGeoTranslation(+16.14,-27.43,   0.0));
  voMagnet -> AddNode(voCoilH,      4, new TGeoTranslation(-16.14,-27.43,   0.0));
  voMagnet -> AddNode(voCoilV,      1, new TGeoTranslation(+16.14,  1.20,+141.8));
  voMagnet -> AddNode(voCoilV,      2, new TGeoTranslation(-16.14,  1.20,+141.8));
  voMagnet -> AddNode(voCoilV,      3, new TGeoTranslation(+16.14,  1.20,-141.8));
  voMagnet -> AddNode(voCoilV,      4, new TGeoTranslation(-16.14,  1.20,-141.8));

  return (TGeoVolume*) voMagnet;
*/
}

//_____________________________________________________________________________
void AliDIPOv2::CreateMaterials()
{
  //
  // Create Materials for Magnetic Dipole version 2
  //
  
  Int_t isxfld1   = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
  Int_t isxfld2   = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->PrecInteg();
  Float_t sxmgmx  = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();
  
  
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






