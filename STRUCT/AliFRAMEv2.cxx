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

//------------------------------------------------------------------------
//  AliFRAMEv2.cxx
//  symmetric space frame with possibility for holes
//  Author: A.Morsch
//------------------------------------------------------------------------

#include <TSystem.h>
#include <TVirtualMC.h>
#include <TString.h>
#include <TGeoManager.h>
#include "AliFRAMEv2.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliConst.h"
#include "AliMC.h"
#include "AliLog.h"
 
ClassImp(AliFRAMEv2)

 
//_____________________________________________________________________________
AliFRAMEv2::AliFRAMEv2():
    fHoles(0)
{
// Constructor
}

//_____________________________________________________________________________
AliFRAMEv2::AliFRAMEv2(const char *name, const char *title)
    : AliFRAME(name,title), 
      fHoles(0)
{
// Constructor
}

//___________________________________________
void AliFRAMEv2::CreateGeometry()
{
//Begin_Html
/*
<img src="picts/frame.gif">
*/
//End_Html


//Begin_Html
/*
<img src="picts/tree_frame.gif">
*/
//End_Html

  Int_t idrotm[2299];
  Int_t *idtmed = fIdtmed->GetArray()-1999;
//
// The Space frame
//
//
  Float_t pbox[3], ptrap[11], ptrd1[4], ppgon[10];
  
  Float_t dx, dy, dz;
  Int_t i, j, jmod;
  jmod = 0;
//
// Constants 
  const Float_t kEps   = 0.01;
  const Int_t kAir   = idtmed[2004];
  const Int_t kSteel = idtmed[2064];
  
  const Float_t krad2deg = 180./TMath::Pi();
  const Float_t kdeg2rad = 1./krad2deg;

  Float_t iFrH   = 114.40;
  Float_t ringH  =   4.00;
  Float_t ringW  =  10.00;
  Float_t longH  =   5.39;
  Float_t longW  =   6.00;  
  Float_t dwl    =   3.14;
  Float_t dwh    =   0.96;

// 
  Float_t dymodU[3] = {70.0, 224.0, 341.};
//  new ?
  Float_t dymodL[3] = {54.0, 178.5, 341.};

//
// Frame mother volume
//
//  ptube[0] = 280.;
//  ptube[1] = 430.;
//  ptube[2] = 376.;
  ppgon[0] =   0.;
  ppgon[1] = 360.;
  ppgon[2] =  18.;
  
  ppgon[3] =   2.;
  ppgon[4] = -376.;
  ppgon[5] =  280.;
  ppgon[6] =  421.;
  
  ppgon[7] =  -ppgon[4]; 
  ppgon[8] =   ppgon[5];
  ppgon[9] =   ppgon[6];

  gMC->Gsvolu("B077", "PGON", kAir, ppgon, 10);
  gMC->Gspos("B077", 1, "ALIC", 0., 0., 0., 0, "ONLY");
//
// Reference plane for TRD
//
  ppgon[6] = ppgon[5] + 0.1;
  ppgon[9] = ppgon[6];
  gMC->Gsvolu("BREF", "PGON", kAir, ppgon, 10);
  gMC->Gspos("BREF", 1, "B077", 0., 0., 0., 0, "ONLY");
  gGeoManager->GetVolume("BREF")->SetVisibility(kFALSE);
//
//  The outer Frame
//

  Float_t dol = 8.75;
  Float_t doh = 5.;
  Float_t ds  = 0.35;
//
// Mother volume
//
  ppgon[0] =   0.;
  ppgon[1] = 360.;
  ppgon[2] =  18.;

  ppgon[3] =   2.;

  ppgon[4] = -350.;
  ppgon[5] =  399.;
  ppgon[6] =  420.7122;
  
  ppgon[7] =  -ppgon[4]; 
  ppgon[8] =   ppgon[5];
  ppgon[9] =   ppgon[6];
  gMC->Gsvolu("B076", "PGON", kAir, ppgon, 10);
  gMC->Gspos("B076", 1, "B077", 0., 0., 0., 0, "ONLY");
//  
// Rings    
//
  dz = 2.*410.2*TMath::Sin(10.*kdeg2rad)-2.*dol*TMath::Cos(10.*kdeg2rad)-
       2.*doh*TMath::Tan(10.*kdeg2rad);
  Float_t l1 = dz/2.;
  Float_t l2 = dz/2.+2.*doh*TMath::Tan(10.*kdeg2rad);

  ptrd1[0] =  l1;
  ptrd1[1] =  l2;
  
  ptrd1[2] =  dol;
  ptrd1[3] =  doh;  
  gMC->Gsvolu("B042", "TRD1", kSteel, ptrd1, 4);

  ptrd1[0] =  ptrd1[0]+ds*(l2-l1)/2./doh;
  ptrd1[1] =  ptrd1[1]-ds*(l2-l1)/2./doh;
  ptrd1[2] =  dol-ds;
  ptrd1[3] =  doh-ds;  
  gMC->Gsvolu("B043", "TRD1", kAir, ptrd1, 4);
  gMC->Gspos("B043", 1, "B042", 0., 0., 0., 0, "ONLY");
//
// longitudinal bars
//
// 170x200x5
//
  pbox[0] = dol;
  pbox[1] = doh;
  pbox[2] = 350.;
  gMC->Gsvolu("B033", "BOX", kSteel, pbox, 3);
  pbox[0] = dol-ds;
  pbox[1] = doh-ds;
  gMC->Gsvolu("B034", "BOX", kAir, pbox, 3);
  gMC->Gspos("B034", 1, "B033", 0., 0., 0., 0, "ONLY");

  pbox[0] =   1.0;
  pbox[1] =   5.0;
  pbox[2] = 375.5;
  

  gMC->Gsvolu("B080", "BOX", kSteel, pbox, 3);
  gMC->Gspos("B080", 1, "B077",  286.01, 0., 0., 0, "ONLY");
  gMC->Gspos("B080", 2, "B077", -286.01, 0., 0., 0, "ONLY");

//
// Diagonal bars (1) 
//
  Float_t h, d, dq, x, theta;
  
  h  = (dymodU[1]-dymodU[0]-2.*dol)*.999;
  d  = 2.*dol;
  dq = h*h+dz*dz;

  x  =  TMath::Sqrt((dz*dz-d*d)/dq + d*d*h*h/dq/dq)+d*h/dq;
  

  theta = krad2deg * TMath::ACos(x);
  
  ptrap[0]  = dz/2.;
  ptrap[1]  = theta;
  ptrap[2]  = 0.;
  ptrap[3]  = doh;
  ptrap[4]  = dol/x;
  ptrap[5]  = ptrap[4];
  ptrap[6]  = 0;
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  ptrap[10] = 0;

  gMC->Gsvolu("B047", "TRAP", kSteel, ptrap, 11);
  ptrap[3]  = doh-ds;
  ptrap[4]  = (dol-ds)/x;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  gMC->Gsvolu("B048", "TRAP", kAir, ptrap, 11);
  gMC->Gspos("B048", 1, "B047", 0.0, 0.0, 0., 0, "ONLY");

/*
 Crosses (inner most) 
       \\  //
        \\//
        //\\
       //  \\
*/
  h  = (2.*dymodU[0]-2.*dol)*.999;
// 
// Mother volume
//
  pbox[0] = h/2;
  pbox[1] = doh;
  pbox[2] = dz/2.;
  gMC->Gsvolu("BM49", "BOX ", kAir, pbox, 3);
  
  
  dq = h*h+dz*dz;
  x  =  TMath::Sqrt((dz*dz-d*d)/dq + d*d*h*h/dq/dq)+d*h/dq;
  theta = krad2deg * TMath::ACos(x);

  ptrap[0]  = dz/2.-kEps;
  ptrap[1]  = theta;
  ptrap[2]  = 0.;
  ptrap[3]  = doh-kEps;
  ptrap[4]  = dol/x;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];

  gMC->Gsvolu("B049", "TRAP", kSteel, ptrap, 11);
  ptrap[0]  = ptrap[0]-kEps;
  ptrap[3]  = (doh-ds);
  ptrap[4]  = (dol-ds)/x;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  gMC->Gsvolu("B050", "TRAP", kAir, ptrap, 11);
  gMC->Gspos("B050", 1, "B049", 0.0, 0.0, 0., 0, "ONLY");
  gMC->Gspos("B049", 1, "BM49", 0.0, 0.0, 0., 0, "ONLY");


  Float_t dd1    = d*TMath::Tan(theta*kdeg2rad);
  Float_t dd2    = d/TMath::Tan(2.*theta*kdeg2rad);
  Float_t theta2 = TMath::ATan(TMath::Abs(dd2-dd1)/d/2.);
  

  ptrap[0] = dol;
  ptrap[1] = theta2*krad2deg;
  ptrap[2] = 0.;
  ptrap[3] = doh;
  ptrap[4] = (dz/2./x-dd1-dd2)/2.;
  ptrap[5] = ptrap[4];
  ptrap[6] = 0.;
  ptrap[7] = ptrap[3];
  ptrap[8] = dz/4./x;
  ptrap[9] = ptrap[8];


  gMC->Gsvolu("B051", "TRAP", kSteel, ptrap, 11);
  Float_t ddx0 = ptrap[8];
  
  Float_t dd1s    = dd1*(1.-2.*ds/d);
  Float_t dd2s    = dd2*(1.-2.*ds/d); 
  Float_t theta2s = TMath::ATan(TMath::Abs(dd2s-dd1s)/(d-2.*ds)/2.);


  ptrap[0] = dol-ds;
  ptrap[1] = theta2s*krad2deg;
  ptrap[2] = 0.;
  ptrap[3] = doh-ds;
  ptrap[4] = ptrap[4]+ds/d/2.*(dd1+dd2);
  ptrap[5] = ptrap[4];
  ptrap[6] = 0.;
  ptrap[7] = ptrap[3];
  ptrap[8] = ptrap[8]-ds/2./d*(dd1+dd2);
  ptrap[9] = ptrap[8];
  
  gMC->Gsvolu("B052", "TRAP", kAir, ptrap, 11);
  gMC->Gspos("B052", 1, "B051", 0.0, 0.0, 0., 0, "ONLY");

  Float_t ddx, ddz, drx, drz, rtheta;

  AliMatrix(idrotm[2001], -theta+180, 0.0, 90.0, 90.0, 90.-theta, 0.0);
  rtheta = (90.-theta)*kdeg2rad;
  ddx = -ddx0-dol*TMath::Tan(theta2);
  ddz = -dol;
  
  drx = TMath::Cos(rtheta) * ddx +TMath::Sin(rtheta) *ddz+pbox[0];
  drz = -TMath::Sin(rtheta) * ddx +TMath::Cos(rtheta) *ddz-pbox[2];
  gMC->Gspos("B051", 1, "BM49", 
	     drx, 0.0, drz,
	     idrotm[2001], "ONLY");

  AliMatrix(idrotm[2002], -theta, 0.0, 90.0, 90.0, 270.-theta, 0.0);
  rtheta = (270.-theta)*kdeg2rad;
  
  drx =  TMath::Cos(rtheta) * ddx +  TMath::Sin(rtheta) * ddz-pbox[0];
  drz = -TMath::Sin(rtheta) * ddx +  TMath::Cos(rtheta) * ddz+pbox[2];
  gMC->Gspos("B051", 2, "BM49", 
	     drx, 0.0, drz,
	     idrotm[2002], "ONLY");

//
// Diagonal bars (3) 
//
  h  = ((dymodU[2]-dymodU[1])-2.*dol)*.999;
  dq = h*h+dz*dz;
  x  =  TMath::Sqrt((dz*dz-d*d)/dq + d*d*h*h/dq/dq)+d*h/dq;
  theta = krad2deg * TMath::ACos(x);
  
  ptrap[0]  = dz/2.;
  ptrap[1]  = theta;
  ptrap[3]  =  doh;
  ptrap[4]  =  dol/x;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];

  gMC->Gsvolu("B045", "TRAP", kSteel, ptrap, 11);
  ptrap[3]  =  doh-ds;
  ptrap[4]  =  (dol-ds)/x;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  gMC->Gsvolu("B046", "TRAP", kAir, ptrap, 11);
  gMC->Gspos("B046", 1, "B045", 0.0, 0.0, 0., 0, "ONLY");

//
// Positioning of diagonal bars
//
// Matrices have been imported from Euclid. Some simplification
// seems possible
//
  AliMatrix(idrotm[2003],   0.0, 0.0, 90.0, 130.0, 90.0,  40.0);
  AliMatrix(idrotm[2004], 180.0, 0.0, 90.0, 130.0, 90.0,  40.0);
  AliMatrix(idrotm[2005], 180.0, 0.0, 90.0, 150.0, 90.0, 240.0);
  AliMatrix(idrotm[2006],   0.0, 0.0, 90.0, 150.0, 90.0, 240.0);
  AliMatrix(idrotm[2007],   0.0, 0.0, 90.0, 170.0, 90.0,  80.0);
  AliMatrix(idrotm[2008], 180.0, 0.0, 90.0, 190.0, 90.0, 280.0);
  AliMatrix(idrotm[2009], 180.0, 0.0, 90.0, 170.0, 90.0,  80.0);
  AliMatrix(idrotm[2010],   0.0, 0.0, 90.0, 190.0, 90.0, 280.0);
  AliMatrix(idrotm[2011],   0.0, 0.0, 90.0, 350.0, 90.0, 260.0);
  AliMatrix(idrotm[2012], 180.0, 0.0, 90.0, 350.0, 90.0, 260.0);
  AliMatrix(idrotm[2013], 180.0, 0.0, 90.0,  10.0, 90.0, 100.0);
  AliMatrix(idrotm[2014],   0.0, 0.0, 90.0,  10.0, 90.0, 100.0);
  AliMatrix(idrotm[2015],   0.0, 0.0, 90.0,  30.0, 90.0, 300.0);
  AliMatrix(idrotm[2016], 180.0, 0.0, 90.0,  30.0, 90.0, 300.0);
  AliMatrix(idrotm[2017], 180.0, 0.0, 90.0,  50.0, 90.0, 140.0);
  AliMatrix(idrotm[2018],   0.0, 0.0, 90.0,  50.0, 90.0, 140.0);

  AliMatrix(idrotm[2019], 180.0, 0.0, 90.0, 130.0, 90.0, 220.0);
  AliMatrix(idrotm[2020], 180.0, 0.0, 90.0,  50.0, 90.0, 320.0);
  AliMatrix(idrotm[2021], 180.0, 0.0, 90.0, 150.0, 90.0,  60.0);
  AliMatrix(idrotm[2022], 180.0, 0.0, 90.0,  30.0, 90.0, 120.0);
  AliMatrix(idrotm[2023], 180.0, 0.0, 90.0, 170.0, 90.0, 260.0);
  AliMatrix(idrotm[2024], 180.0, 0.0, 90.0, 190.0, 90.0, 100.0);
  AliMatrix(idrotm[2025], 180.0, 0.0, 90.0, 350.0, 90.0,  80.0);
  AliMatrix(idrotm[2026], 180.0, 0.0, 90.0,  10.0, 90.0, 280.0);
  
  AliMatrix(idrotm[2027],   0.0, 0.0, 90.0,  50.0, 90.0, 320.0);
  AliMatrix(idrotm[2028],   0.0, 0.0, 90.0, 150.0, 90.0,  60.0); 
  AliMatrix(idrotm[2029],   0.0, 0.0, 90.0,  30.0, 90.0, 120.0);
  AliMatrix(idrotm[2030],   0.0, 0.0, 90.0,  10.0, 90.0, 280.0);
  AliMatrix(idrotm[2031],   0.0, 0.0, 90.0, 170.0, 90.0, 260.0);
  AliMatrix(idrotm[2032],   0.0, 0.0, 90.0, 190.0, 90.0, 100.0);
  AliMatrix(idrotm[2033],   0.0, 0.0, 90.0, 350.0, 90.0,  80.0);
  
  Float_t rd =  410.56;
  dz = (dymodU[1]+dymodU[0])/2.;
  Float_t dz2 =  (dymodU[1]+dymodU[2])/2.;

//
//  phi = 40
//
  Float_t  phi = 40;
  dx = rd * TMath::Sin(phi*kdeg2rad);
  dy = rd * TMath::Cos(phi*kdeg2rad);
  

  gMC->Gspos("B045", 1, "B076", -dx,  dy,  dz2, idrotm[2019], "ONLY");
  gMC->Gspos("B045", 2, "B076", -dx,  dy, -dz2, idrotm[2003], "ONLY"); // ?
  gMC->Gspos("B045", 3, "B076",  dx,  dy,  dz2, idrotm[2020], "ONLY");
  gMC->Gspos("B045", 4, "B076",  dx,  dy, -dz2, idrotm[2027], "ONLY");


//
//  phi = 60
//

  phi = 60;
  dx = rd * TMath::Sin(phi*kdeg2rad);
  dy = rd * TMath::Cos(phi*kdeg2rad);

  gMC->Gspos("B045", 5, "B076", -dx,  dy,  dz2, idrotm[2021], "ONLY");
  gMC->Gspos("B045", 6, "B076", -dx,  dy, -dz2, idrotm[2028], "ONLY");
  gMC->Gspos("B045", 7, "B076",  dx,  dy,  dz2, idrotm[2022], "ONLY");
  gMC->Gspos("B045", 8, "B076",  dx,  dy, -dz2, idrotm[2029], "ONLY");

//
//  phi = 80
//

  phi = 80;
  dx = rd * TMath::Sin(phi*kdeg2rad);
  dy = rd * TMath::Cos(phi*kdeg2rad);

  gMC->Gspos("B047", 13, "B076", -dx, -dy,  dz, idrotm[2008], "ONLY");
  gMC->Gspos("B047", 14, "B076", -dx, -dy, -dz, idrotm[2010], "ONLY");
  gMC->Gspos("B047", 15, "B076",  dx, -dy,  dz, idrotm[2012], "ONLY");
  gMC->Gspos("B047", 16, "B076",  dx, -dy, -dz, idrotm[2011], "ONLY");

  gMC->Gspos("B045",  9, "B076", -dx,  dy,  dz2, idrotm[2023], "ONLY");
  gMC->Gspos("B045", 10, "B076", -dx,  dy, -dz2, idrotm[2031], "ONLY");
  gMC->Gspos("B045", 11, "B076",  dx,  dy,  dz2, idrotm[2026], "ONLY");
  gMC->Gspos("B045", 12, "B076",  dx,  dy, -dz2, idrotm[2030], "ONLY");

  gMC->Gspos("B045", 13, "B076", -dx, -dy,  dz2, idrotm[2024], "ONLY");
  gMC->Gspos("B045", 14, "B076", -dx, -dy, -dz2, idrotm[2032], "ONLY");
  gMC->Gspos("B045", 15, "B076",  dx, -dy,  dz2, idrotm[2025], "ONLY");
  gMC->Gspos("B045", 16, "B076",  dx, -dy, -dz2, idrotm[2033], "ONLY");

  gMC->Gspos("BM49", 7, "B076",  dx, -dy,  0., idrotm[2025], "ONLY");
  gMC->Gspos("BM49", 8, "B076", -dx, -dy,  0., idrotm[2024], "ONLY");


// The internal frame
//
//
//
//  Mother Volumes
//
  ptrd1[0] = 50.18;
  ptrd1[1] = 70.35;
  ptrd1[2] = 375.5;
  ptrd1[3] =  57.2;  
  
  Float_t r      = 341.8;
  Float_t rout1  = 410.564;
  Float_t rout2  = 415.2;
  TString module[18];
  
  for (i = 0; i < 18; i++) {
      // Create volume i 
      char name[16];
      Int_t mod = i + 13;
      if (mod > 17) mod -= 18;
      sprintf(name, "BSEGMO%d", mod);
      gMC->Gsvolu(name, "TRD1", kAir, ptrd1, 4);
      gGeoManager->GetVolume(name)->SetVisibility(kFALSE);

      module[i] = name;
      // Place volume i
      Float_t phi  = i * 20.;
      Float_t phi2 = 270 + phi;
      if (phi2 >= 360.) phi2 -= 360.;
      
      dx =  TMath::Sin(phi*kdeg2rad)*r;
      dy = -TMath::Cos(phi*kdeg2rad)*r;
      
      AliMatrix(idrotm[2034+i],  90.0, phi, 0., 0., 90., phi2);  
      gMC->Gspos(name, 1, "B077", dx, dy, 0., idrotm[2034+i], "ONLY");
//
//    Position elements of outer Frame
//
      dx =  TMath::Sin(phi*kdeg2rad)*rout1;
      dy = -TMath::Cos(phi*kdeg2rad)*rout1;
      for (j = 0; j < 3; j++)
      {
	  dz = dymodU[j];
	  gMC->Gspos("B042", 6*i+2*j+1, "B076", dx, dy,  dz, idrotm[2034+i], "ONLY");	  
	  gMC->Gspos("B042", 6*i+2*j+2, "B076", dx, dy, -dz, idrotm[2034+i], "ONLY");	  
      }

      phi = i*20.+10;
      phi2 = 270+phi;
      AliMatrix(idrotm[2052+i],  90.0, phi, 90., phi2, 0., 0.);  

      dx =  TMath::Sin(phi*kdeg2rad)*rout2;
      dy = -TMath::Cos(phi*kdeg2rad)*rout2;
      gMC->Gspos("B033", i+1, "B076", dx, dy,  0., idrotm[2052+i], "ONLY");	  
//
  }
// Internal Frame rings
//
//
// new specs: 40x100x6 for inner rings
//            30x135x6 for front and rear rings
//
// currently no distinction between front/rear and inner rings
// 
//
//
  pbox[0] = 50.;
  pbox[1] =  ringW/2.;
  pbox[2] =  ringH/2.;
  
  gMC->Gsvolu("B072", "BOX ", kSteel, pbox, 3);

  pbox[1] =  pbox[1] - 0.6;
  pbox[2] =  pbox[2] - 0.6;  
  gMC->Gsvolu("B073", "BOX ", kAir, pbox, 3);
  gMC->Gspos("B073", 1, "B072", 0., 0., 0., 0, "ONLY");

// Web frame 0-degree
//
// h x w x s = 60x40x4 
// (attention: element is are half bars, "U" shaped)  
//
  dz  = 106.2;
  d   = 2.*dwl;
  h   = dymodU[0]-dymodL[0];
  dq  = h*h+dz*dz;

  x  =  TMath::Sqrt((dz*dz-d*d)/dq + d*d*h*h/dq/dq)+d*h/dq;
  
  theta = TMath::ACos(x);
  Float_t thetaP = (TMath::Pi()-theta)*krad2deg;
  
  ptrap[0]  = dz/2.;
  ptrap[1]  = thetaP;
  ptrap[2]  =  0.;
  ptrap[3]  =  dwh;
  ptrap[4]  =  dwl/cos(theta);
  ptrap[5]  = ptrap[4];
  ptrap[6]  = 0;
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  ptrap[10] = 0;
  gMC->Gsvolu("B063", "TRAP", kSteel, ptrap, 11);
  ptrap[3]  =  dwh - 0.2;
  ptrap[4]  =  (dwl-0.4)/cos(theta);
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  gMC->Gsvolu("B064", "TRAP", kAir, ptrap, 11);
  gMC->Gspos("B064", 1, "B063", 0.0, -0.2, 0., 0, "ONLY");


  h  = 106.2;
  d  = 2.*dwl;
  dz = dymodU[1]-dymodL[1];
  dq = h*h+dz*dz;

  x  =  TMath::Sqrt((dz*dz-d*d)/dq + d*d*h*h/dq/dq)+d*h/dq;
  

  theta = krad2deg * TMath::ACos(x);
  
  ptrap[0]  = dz/2.;
  ptrap[1]  = theta;
  ptrap[2]  =  0.;
  ptrap[3]  =  dwh;
  ptrap[4]  =  dwl/x;
  ptrap[5]  = ptrap[4];
  ptrap[6]  = 0;
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  ptrap[10] = 0;
  gMC->Gsvolu("B163", "TRAP", kSteel, ptrap, 11);
  ptrap[3]  =  dwh - 0.2;
  ptrap[4]  =  (dwl-0.4)/x;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  gMC->Gsvolu("B164", "TRAP", kAir, ptrap, 11);
  gMC->Gspos("B164", 1, "B163", 0.0, -0.2, 0., 0, "ONLY");
 
 
  pbox[0] = dwh;
  pbox[1] = dwl;
  pbox[2] = (iFrH-ringH-longH)/2.;
  gMC->Gsvolu("B263", "BOX ", kSteel, pbox, 3);
  pbox[0] = dwh-0.2;
  pbox[1] = dwl-0.4;
  
  gMC->Gsvolu("B264", "BOX ", kAir, pbox, 3);
  gMC->Gspos("B264", 1, "B263", 0.2, 0., 0., 0, "ONLY");

  
  AliMatrix(idrotm[2070],  90.0,   0.0,  90.0, 270.0,   0.0,   0.0);  
//
  AliMatrix(idrotm[2071], 100.0,   0.0,  90.0, 270.0,  10.0,   0.0);
  AliMatrix(idrotm[2072], 100.0,   0.0,  90.0,  90.0,  10.0,   0.0);
  AliMatrix(idrotm[2073], 100.0, 180.0,  90.0, 270.0,  10.0, 180.0);
  AliMatrix(idrotm[2074], 100.0, 180.0,  90.0,  90.0,  10.0, 180.0);
//
  AliMatrix(idrotm[2075],  10.0,   0.0,  80.0, 180.0,  90.0,  90.0);
  AliMatrix(idrotm[2076],  10.0,   0.0,  80.0, 180.0,  90.0, 270.0);
  AliMatrix(idrotm[2077],  10.0, 180.0,  80.0,   0.0,  90.0,  90.0);
  AliMatrix(idrotm[2078],  10.0, 180.0,  80.0,   0.0,  90.0, 270.0);
//
  AliMatrix(idrotm[2079], 170.0, 180.0,  80.0, 180.0,  90.0,  90.0);
  AliMatrix(idrotm[2080], 170.0, 180.0,  80.0, 180.0,  90.0, 270.0);
  AliMatrix(idrotm[2081], 170.0,   0.0,  80.0,   0.0,  90.0,  90.0);
  AliMatrix(idrotm[2082], 170.0,   0.0,  80.0,   0.0,  90.0, 270.0);

//
  AliMatrix(idrotm[2083], 170.0,   0.0,  90.0,  90.0,  80.0,   0.0);
  AliMatrix(idrotm[2084], 170.0, 180.0,  90.0,  90.0,  80.0, 180.0);
  AliMatrix(idrotm[2085],  90.0, 180.0,  90.0,  90.0,   0.0,   0.0);
//
  AliMatrix(idrotm[2086],  90.0,  90.0, 100.0,    0., 170.0, 180.0);
  AliMatrix(idrotm[2087],  90.0,  90.0, 100.0,    0.,  10.0,   0.0);
  AliMatrix(idrotm[2088],  90.0,  90.0,  80.0,    0., 170.0,   0.0);
  AliMatrix(idrotm[2089],  90.0,  90.0,  80.0,    0., -10.0,   0.0);
//
  AliMatrix(idrotm[2090],  90.0,   0.0,   0.0,    0.,   90.0, 90.0);
  AliMatrix(idrotm[2091],   0.0,   0.0,  90.0,   90.,   90.0,  0.0);
//
// web frame diagonal (outer)
//  
  h  = 106.2;
  d  = 2.*dwl;
  dz = dymodL[2]-dymodU[1]-dwl;
  dq = h*h+dz*dz;

  x  =  TMath::Sqrt((dz*dz-d*d)/dq + d*d*h*h/dq/dq)+d*h/dq;
  

  theta = krad2deg * TMath::ACos(x);
  
  ptrap[0]  = dz/2.;
  ptrap[1]  = theta;
  ptrap[2]  =  0.;
  ptrap[3]  =  dwh;
  ptrap[4]  =  dwl/x;
  ptrap[5]  = ptrap[4];
  ptrap[6]  = 0;
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  ptrap[10] = 0;
  gMC->Gsvolu("B065", "TRAP", kSteel, ptrap, 11);
  ptrap[3]  =  dwh - 0.2;
  ptrap[4]  =  (dwl-0.4)/x;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  gMC->Gsvolu("B066", "TRAP", kAir, ptrap, 11);
  gMC->Gspos("B066", 1, "B065", 0.0, -0.2, 0., 0, "ONLY");

//
// web frame diagonal (inner)
//
  dz = dymodL[1]-dymodU[0];
  dq = h*h+dz*dz;
  x  =  TMath::Sqrt((dz*dz-d*d)/dq + d*d*h*h/dq/dq)+d*h/dq;
  

  theta = krad2deg * TMath::ACos(x);
  
  ptrap[0]  = (dymodL[1]-dymodU[0])/2.;
  ptrap[1]  = theta;
  ptrap[2]  =  0.;
  ptrap[3]  =  dwh;
  ptrap[4]  =  dwl/x;
  ptrap[5]  = ptrap[4];
  ptrap[6]  = 0;
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  ptrap[10] = 0;
  gMC->Gsvolu("B067", "TRAP", kSteel, ptrap, 11);
  ptrap[3]  =  dwh-0.2;
  ptrap[4]  =  (dwl-0.4)/x;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  gMC->Gsvolu("B068", "TRAP", kAir, ptrap, 11);
  gMC->Gspos("B068", 1, "B067", 0.0, -0.2, 0., 0, "ONLY");


  dz = -iFrH/2.+ringH/2.+kEps;
  
  for (jmod = 0; jmod< 18; jmod++)
  {
// ring bars
      for (i = 0; i < 3; i++) {
	  gMC->Gspos("B072", 6*jmod+i+1, module[jmod], 0,  dymodL[i], dz, 0, "ONLY");
	  gMC->Gspos("B072", 6*jmod+i+4, module[jmod], 0, -dymodL[i], dz, idrotm[2070], "ONLY");      

// 0-deg web
	  if (i == 2) {
	      gMC->Gspos("B263", 4*jmod+1,  module[jmod],  60.0732,  dymodU[2], 4.6669, idrotm[2072], "ONLY");
	      gMC->Gspos("B263", 4*jmod+2,  module[jmod],  60.0732, -dymodU[2], 4.6669, idrotm[2071], "ONLY");      
	      gMC->Gspos("B263", 4*jmod+3,  module[jmod], -60.0732,  dymodU[2], 4.6669, idrotm[2074], "ONLY");
	      gMC->Gspos("B263", 4*jmod+4,  module[jmod], -60.0732, -dymodU[2], 4.6669, idrotm[2073], "ONLY");      
	  }
      }
  }
  
// outer diagonal web

  dy = (dymodU[0]+dymodL[0])/2.;
  for (jmod = 0; jmod < 18; jmod++) {
      gMC->Gspos("B063", 4*jmod+1, module[jmod],  60.0732,   dy, 4.6669, idrotm[2086], "ONLY");
      gMC->Gspos("B063", 4*jmod+2, module[jmod],  60.0732,  -dy, 4.6669, idrotm[2087], "ONLY");
      gMC->Gspos("B063", 4*jmod+3, module[jmod], -60.0732,   dy, 4.6669, idrotm[2088], "ONLY");
      gMC->Gspos("B063", 4*jmod+4, module[jmod], -60.0732,  -dy, 4.6669, idrotm[2089], "ONLY");
  }

  dy = (dymodU[1]+dymodL[1])/2.;
  for (jmod = 0; jmod < 18; jmod++) {
      gMC->Gspos("B163", 4*jmod+1, module[jmod],  60.0732,   dy, 4.6669, idrotm[2080], "ONLY");
      gMC->Gspos("B163", 4*jmod+2, module[jmod],  60.0732,  -dy, 4.6669, idrotm[2079], "ONLY");
      gMC->Gspos("B163", 4*jmod+3, module[jmod], -60.0732,   dy, 4.6669, idrotm[2082], "ONLY");
      gMC->Gspos("B163", 4*jmod+4, module[jmod], -60.0732,  -dy, 4.6669, idrotm[2081], "ONLY");
  }


  dy = (dymodL[2]+dymodU[1])/2.-dwl/2.;
  for (jmod = 0; jmod < 18; jmod++) {
      gMC->Gspos("B065", 4*jmod+1, module[jmod],  60.0732,   dy, 4.6669, idrotm[2076], "ONLY");
      gMC->Gspos("B065", 4*jmod+2, module[jmod],  60.0732,  -dy, 4.6669, idrotm[2075], "ONLY");
      gMC->Gspos("B065", 4*jmod+3, module[jmod], -60.0732,   dy, 4.6669, idrotm[2078], "ONLY");
      gMC->Gspos("B065", 4*jmod+4, module[jmod], -60.0732,  -dy, 4.6669, idrotm[2077], "ONLY");
  }
  

  dy = (dymodL[1]+dymodU[0])/2.;

  for (jmod = 0; jmod < 18; jmod++) {
      gMC->Gspos("B067", 4*jmod+1, module[jmod],  60.0732,   dy, 4.6669, idrotm[2076], "ONLY");
      gMC->Gspos("B067", 4*jmod+2, module[jmod],  60.0732,  -dy, 4.6669, idrotm[2075], "ONLY");
      gMC->Gspos("B067", 4*jmod+3, module[jmod], -60.0732,   dy, 4.6669, idrotm[2078], "ONLY");
      gMC->Gspos("B067", 4*jmod+4, module[jmod], -60.0732,  -dy, 4.6669, idrotm[2077], "ONLY");
  }
 
// longitudinal bars (TPC rails attached)
//  new specs:
//  h x w x s = 100 x 75 x 6 
//  current: 
//  Attention: 2 "U" shaped half rods per cell 
//
//  not yet used 
//
  ptrap[0]  =   2.50;
  ptrap[1]  =  10.00;
  ptrap[2]  =   0.00;
  ptrap[3]  = 350.00;
  ptrap[4]  =   3.75;
  ptrap[5]  = ptrap[4];
  ptrap[6]  = 0;
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  ptrap[10] = 0;
//  gMC->Gsvolu("B059", "TRAP", kSteel, ptrap, 11);
  ptrap[0]  =  2.2;
  ptrap[4]  =  2.15;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  //gMC->Gsvolu("B062", "TRAP", kAir, ptrap, 11);
  //gMC->Gspos("B062", 1, "B059", 0.0, 0., 0., 0, "ONLY");

//
// longitudinal bars (no TPC rails attached)
// new specs: h x w x s = 60 x 60 x 3
// (was: 75 x 100 x 5?)
//
//
// 
  ptrap[0]  = longW/4.;
  ptrap[4]  = longH/2.;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];

  gMC->Gsvolu("BA59", "TRAP", kSteel, ptrap, 11);
  ptrap[0]  = longW/4.-0.15;
  ptrap[4]  = longH/2.-0.30;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  gMC->Gsvolu("BA62", "TRAP", kAir, ptrap, 11);
  gMC->Gspos("BA62", 1, "BA59", 0.0, 0.0, -0.15, 0, "ONLY");

  dz = -iFrH/2.+ringH+longH/2.;

  for (jmod = 0; jmod < 18; jmod++) {
      gMC->Gspos("BA59", 2*jmod+1, module[jmod],  49.6476, 0.0, dz, idrotm[2084], "ONLY");
      gMC->Gspos("BA59", 2*jmod+2, module[jmod], -49.6476, 0.0, dz, idrotm[2083], "ONLY");
  }

  
//
// Thermal shield
//
  Float_t dyM  =  99.0;
  MakeHeatScreen("M",   dyM, idrotm[2090], idrotm[2091]);
  Float_t dyAM = 119.5;
  MakeHeatScreen("AM", dyAM, idrotm[2090], idrotm[2091]);
  Float_t dyA  = 128.0;
  MakeHeatScreen("A" ,  dyA, idrotm[2090], idrotm[2091]);
//
//
//
  dz = -57.2 + 5.2 - 0.7;  
  for (i = 0; i < 18; i++) {

      char nameMo[16];
      sprintf(nameMo, "BSEGMO%d",i);
      // M
      gMC->Gspos("BTSH_M" , i+1 , nameMo,  0., 0., dz, 0, "ONLY"); 
      // AM, CM
      dy = dymodL[0] + dyAM / 2.;
      gMC->Gspos("BTSH_AM", i+ 1, nameMo, 0.,  dy, dz, 0, "ONLY"); 
      gMC->Gspos("BTSH_AM", i+19, nameMo, 0., -dy, dz, 0, "ONLY"); 
      // A, C
      dy = dymodL[1] + dyA / 2.;
      gMC->Gspos("BTSH_A" , i+ 1, nameMo, 0.,  dy, dz, 0, "ONLY"); 
      gMC->Gspos("BTSH_A" , i+19, nameMo, 0., -dy, dz, 0, "ONLY"); 
}
  
  
  //
  // TRD mother volumes
  //

  ptrd1[0] = 47.4405;   // CBL 28/6/2006
  ptrd1[1] = 61.1765;   // CBL
  ptrd1[2] = 375.5;     // CBL
  ptrd1[3] = 38.95;     // CBL



  for (i = 0; i < 18; i++) {
    char nameCh[16];
    sprintf(nameCh, "BTRD%d",i);
    char nameMo[16];
    sprintf(nameMo, "BSEGMO%d",i);
    gMC->Gsvolu(nameCh, "TRD1", kAir, ptrd1, 4);
    gGeoManager->GetVolume(nameCh)->SetVisibility(kFALSE);
    gMC->Gspos(nameCh, 1, nameMo, 0., 0., -11.75, 0, "ONLY"); // CBL 28/6/2006
  }

//
// TOF  mother volumes
//
  ptrd1[0] = 63.2061;
  ptrd1[1] = 68.3192;
  ptrd1[2] = 375.5;
  ptrd1[3] = 14.5;
  for (i = 0; i < 18; i++) {
    char nameCh[16];
    sprintf(nameCh, "BTOF%d",i);
    char nameMo[16];
    sprintf(nameMo, "BSEGMO%d",i);
    gMC->Gsvolu(nameCh, "TRD1", kAir, ptrd1, 4);
    gGeoManager->GetVolume(nameCh)->SetVisibility(kFALSE);
    gMC->Gspos(nameCh, 1, nameMo, 0., 0., 42.69, 0, "ONLY");
  }

//
//    Geometry of Rails starts here
//
//
//
//    Rails for space-frame
//
  Float_t rbox[3];

  rbox[0] =  25.00;
  rbox[1] =  27.50;
  rbox[2] = 600.00;
  gMC->Gsvolu("BRS1", "BOX", kAir, rbox, 3);
  
  rbox[0] =  25.00;
  rbox[1] =   3.75;
  gMC->Gsvolu("BRS2", "BOX", kSteel, rbox, 3);
  
  rbox[0] =   3.00;
  rbox[1] =  20.00;
  gMC->Gsvolu("BRS3", "BOX", kSteel, rbox, 3);
  
  gMC->Gspos("BRS2", 1, "BRS1", 0., -27.5+3.75, 0., 0, "ONLY");
  gMC->Gspos("BRS2", 2, "BRS1", 0.,  27.5-3.75, 0., 0, "ONLY");
  gMC->Gspos("BRS3", 1, "BRS1", 0.,         0., 0., 0, "ONLY");
  gMC->Gspos("BRS1", 1, "ALIC", -430.-3.,    -190., 0., 0, "ONLY");
  gMC->Gspos("BRS1", 2, "ALIC",  430.+3.,    -190., 0., 0, "ONLY");

  rbox[0] =    3.0;
  rbox[1] =  145./4.;
  rbox[2] =   25.0;
  gMC->Gsvolu("BRS4", "BOX", kSteel, rbox, 3);

  gMC->Gspos("BRS4", 1, "ALIC",  430.+3.,    -190.+55./2.+rbox[1],  224., 0, "ONLY");
  gMC->Gspos("BRS4", 2, "ALIC",  430.+3.,    -190.+55./2.+rbox[1], -224., 0, "ONLY");
//  gMC->Gspos("BRS4", 3, "ALIC", -430.+3,    -180.+55./2.+rbox[1],  224., 0, "ONLY");
//  gMC->Gspos("BRS4", 4, "ALIC", -430.+3,    -180.+55./2.+rbox[1], -224., 0, "ONLY");



  //
  // The Backframe
  //
  // Inner radius 
  Float_t kBFMRin = 270.0;
  // Outer Radius
  Float_t kBFMRou = 417.5;
  // Width
  Float_t kBFMdz  = 118.0;
  //
  //
  // Rings
  Float_t kBFRdr   =  7.5;
  Float_t kBFRdz   =  8.0;
  //
  //
  // Bars and Spokes
  //
  Float_t kBFBd   =   8.0;
  Float_t kBFBdd  =   0.6;
  

  // The Mother volume
  Float_t tpar[3];
  tpar[0] = kBFMRin;
  tpar[1] = kBFMRou;
  tpar[2] = kBFMdz / 2.;
  gMC->Gsvolu("BFMO", "TUBE", kAir, tpar, 3);  
  
  // Rings
  //
  // Inner Ring
  tpar[0] =  kBFMRin;
  tpar[1] =  tpar[0] +  kBFRdr;
  tpar[2] =  kBFRdz / 2.;
  
  gMC->Gsvolu("BFIR", "TUBE", kSteel, tpar, 3);  
  
  tpar[0] =  tpar[0] +  kBFBdd;
  tpar[1] =  tpar[1] -  kBFBdd;
  tpar[2] =  (kBFRdz - 2. * kBFBdd) / 2.;

  gMC->Gsvolu("BFII", "TUBE", kAir, tpar, 3);  
  gMC->Gspos("BFII", 1, "BFIR", 0., 0., 0., 0, "ONLY");  

  //
  // Outer RING
  tpar[0] =  kBFMRou - kBFRdr;
  tpar[1] =  kBFMRou;
  tpar[2] =  kBFRdz / 2.;
  
  gMC->Gsvolu("BFOR", "TUBE", kSteel, tpar, 3);  
  
  tpar[0] =  tpar[0] +  kBFBdd;
  tpar[1] =  tpar[1] -  kBFBdd;
  tpar[2] =  (kBFRdz - 2. * kBFBdd) / 2.;

  gMC->Gsvolu("BFOO", "TUBE", kAir, tpar, 3);  
  gMC->Gspos("BFOO", 1, "BFOR", 0., 0., 0., 0, "ONLY");  


  dz = kBFMdz/2. -  kBFRdz / 2.;
  gMC->Gspos("BFIR", 1, "BFMO", 0., 0.,  dz, 0, "ONLY");  
  gMC->Gspos("BFIR", 2, "BFMO", 0., 0., -dz, 0, "ONLY");  
  gMC->Gspos("BFOR", 1, "BFMO", 0., 0.,  dz, 0, "ONLY");  
  gMC->Gspos("BFOR", 2, "BFMO", 0., 0., -dz, 0, "ONLY");  
  
  // 
  // Longitudinal Bars
  // 
  Float_t bpar[3];
  
  bpar[0] =  kBFBd/2;
  bpar[1] =  bpar[0];
  bpar[2] =  kBFMdz/2.  - kBFBd;
  gMC->Gsvolu("BFLB", "BOX ", kSteel, bpar, 3); 

  bpar[0] = bpar[0] - kBFBdd;
  bpar[1] = bpar[1] - kBFBdd;
  bpar[2] = bpar[2] - kBFBdd;
  gMC->Gsvolu("BFLL", "BOX ", kAir, bpar, 3); 
  gMC->Gspos("BFLL", 1, "BFLB", 0., 0., 0., 0, "ONLY");  

  for (i = 0; i < 18; i++)
  {
      Float_t ro   = kBFMRou - kBFBd / 2.;
      Float_t ri   = kBFMRin + kBFBd / 2.;

      Float_t phi = Float_t(i) * 20.;
      
      Float_t x = ri * TMath::Cos(phi * kDegrad);
      Float_t y = ri * TMath::Sin(phi * kDegrad);
      AliMatrix(idrotm[2090+i],  90.0, phi,  90.0, phi + 270., 0., 0.);
      
      gMC->Gspos("BFLB", i + 1, "BFMO", x, y, 0., idrotm[2090 + i], "ONLY");      

      x = ro * TMath::Cos(phi * kDegrad);
      y = ro * TMath::Sin(phi * kDegrad);

      gMC->Gspos("BFLB", i + 19, "BFMO", x, y, 0., idrotm[2090 +i], "ONLY");       
 }

  // 
  // Radial Bars
  // 
  bpar[0] =  (kBFMRou - kBFMRin - 2. * kBFRdr) / 2.;
  bpar[1] =  kBFBd/2;
  bpar[2] =  bpar[1];
  //
  // Avoid overlap with circle
  Float_t rr    = kBFMRou - kBFRdr;
  Float_t delta = rr - TMath::Sqrt(rr * rr - kBFBd * kBFBd / 4.) + 0.01;
  bpar[0] -= delta /2.;
  

  gMC->Gsvolu("BFRB", "BOX ", kSteel, bpar, 3); 

  bpar[0] = bpar[0] - kBFBdd;
  bpar[1] = bpar[1] - kBFBdd;
  bpar[2] = bpar[2] - kBFBdd;
  gMC->Gsvolu("BFRR", "BOX ", kAir, bpar, 3); 
  gMC->Gspos("BFRR", 1, "BFRB", 0., 0., 0., 0, "ONLY");  

  Int_t iphi[10] = {0, 1, 3, 6, 8, 9, 10, 12, 15, 17};
  
  for (i = 0; i < 10; i++)
  {
      
      Float_t r   = (kBFMRin + kBFMRou)/2.;
      Float_t phi = Float_t(iphi[i]) * 20.;
      
      Float_t x = r * TMath::Cos(phi * kDegrad);
      Float_t y = r * TMath::Sin(phi * kDegrad);
      
      gMC->Gspos("BFRB", i + 1,  "BFMO", x, y,  dz, idrotm[2034 + iphi[i]], "ONLY");      
      gMC->Gspos("BFRB", i + 11, "BFMO", x, y, -dz, idrotm[2034 + iphi[i]], "ONLY");      

 }

  gMC->Gspos("BFMO", i + 19, "ALIC", 0, 0, - 376. - kBFMdz/2. - 0.5 , 0, "ONLY");       



//
//
//  The Baby Frame
//
//
  //
  // Inner radius 
  Float_t kBBMRin = 278.0;
  // Outer Radius
  Float_t kBBMRou = 410.5;
  // Width
  Float_t kBBMdz  = 223.0;
  Float_t kBBBdz  = 6.0;
  Float_t kBBBdd  = 0.6;

  
  // The Mother volume

  ppgon[0] =   0.;
  ppgon[1] = 360.;
  ppgon[2] =  18.;
  
  ppgon[3] =   2.;
  ppgon[4] = -kBBMdz / 2. ;
  ppgon[5] =  kBBMRin;
  ppgon[6] =  kBBMRou;
  
  ppgon[7] =  -ppgon[4]; 
  ppgon[8] =   ppgon[5];
  ppgon[9] =   ppgon[6];

  gMC->Gsvolu("BBMO", "PGON", kAir, ppgon, 10);
  gMC->Gsdvn("BBCE", "BBMO", 18, 2);

  // Longitudinal bars
  bpar[0] =  kBBBdz/2.;
  bpar[1] =  bpar[0];
  bpar[2] =  kBBMdz/2.  - kBBBdz;
  gMC->Gsvolu("BBLB", "BOX ", kSteel, bpar, 3); 
  bpar[0] -= kBBBdd;
  bpar[1] -= kBBBdd;
  bpar[2] -= kBBBdd;
  gMC->Gsvolu("BBLL", "BOX ", kAir, bpar, 3); 
  gMC->Gspos("BBLL", 1, "BBLB", 0., 0., 0., 0, "ONLY"); 

  dx = kBBMRin + kBBBdz/2. + (bpar[1] + kBBBdd) * TMath::Sin(10. * kDegrad);
  dy = dx * TMath::Tan(10. * kDegrad) - kBBBdz/2./TMath::Cos(10. * kDegrad);
  gMC->Gspos("BBLB", 1, "BBCE", dx, dy, 0., idrotm[2052], "ONLY"); 

  dx = kBBMRou - kBBBdz/2. - (bpar[1] + kBBBdd) * TMath::Sin(10. * kDegrad);
  dy = dx * TMath::Tan(10. * kDegrad) - kBBBdz/2./TMath::Cos(10. * kDegrad);
 
  gMC->Gspos("BBLB", 2, "BBCE", dx, dy, 0., idrotm[2052], "ONLY");  

  // 
  // Radial Bars
  // 
  bpar[0] =  (kBBMRou - kBBMRin) / 2. - kBBBdz;
  bpar[1] =  kBBBdz/2;
  bpar[2] =  bpar[1];

  gMC->Gsvolu("BBRB", "BOX ", kSteel, bpar, 3); 
  bpar[0] -= kBBBdd;
  bpar[1] -= kBBBdd;
  bpar[2] -= kBBBdd;
  gMC->Gsvolu("BBRR", "BOX ", kAir, bpar, 3); 
  gMC->Gspos("BBRR", 1, "BBRB", 0., 0., 0., 0, "ONLY"); 


  dx = (kBBMRou + kBBMRin) / 2.;
  dy = ((kBBMRou + kBBMRin)/ 2) *  TMath::Tan(10 * kDegrad) - kBBBdz / 2./ TMath::Cos(10 * kDegrad);
  dz = kBBMdz/2. -  kBBBdz / 2.;

  gMC->Gspos("BBRB", 1, "BBCE", dx, dy,   dz, idrotm[2052], "ONLY");  
  gMC->Gspos("BBRB", 2, "BBCE", dx, dy, - dz, idrotm[2052], "ONLY");  
  gMC->Gspos("BBRB", 3, "BBCE", dx, dy,   0., idrotm[2052], "ONLY");  
 
 //
 // Circular bars 
 //
 //  Inner
  
  bpar[1] =  kBBMRin * TMath::Sin(10. * kDegrad);
  bpar[0] =  kBBBdz/2;
  bpar[2] =  bpar[0];
  gMC->Gsvolu("BBC1", "BOX ", kSteel, bpar, 3); 
  bpar[0] -= kBBBdd;
  bpar[1] -= kBBBdd;
  bpar[2] -= kBBBdd;
  gMC->Gsvolu("BBC2", "BOX ", kAir, bpar, 3); 
  gMC->Gspos("BBC2", 1, "BBC1", 0., 0., 0., 0, "ONLY"); 
  dx = kBBMRin + kBBBdz/2;
  dy = 0.;
  gMC->Gspos("BBC1", 1, "BBCE", dx, dy,   dz, 0, "ONLY");  
  gMC->Gspos("BBC1", 2, "BBCE", dx, dy,  -dz, 0, "ONLY");  
  //
  // Outer
  bpar[1] =  (kBBMRou - kBBBdz) * TMath::Sin(10. * kDegrad);
  bpar[0] =  kBBBdz/2;
  bpar[2] =  bpar[0];
  gMC->Gsvolu("BBC3", "BOX ", kSteel, bpar, 3); 
  bpar[0] -= kBBBdd;
  bpar[1] -= kBBBdd;
  bpar[2] -= kBBBdd;
  gMC->Gsvolu("BBC4", "BOX ", kAir, bpar, 3); 
  gMC->Gspos("BBC4", 1, "BBC3", 0., 0., 0., 0, "ONLY"); 
  dx = kBBMRou - kBBBdz/2;
  dy = 0.;
  gMC->Gspos("BBC3", 1, "BBCE", dx, dy,   dz, 0, "ONLY");  
  gMC->Gspos("BBC3", 2, "BBCE", dx, dy, - dz, 0, "ONLY");
  //
  // Diagonal Bars
  //
  h  = (kBBMRou - kBBMRin - 2. * kBBBdz);;
  d  = kBBBdz;
  dz = kBBMdz/2. - 1.6 * kBBBdz;
  dq = h*h+dz*dz;

  x  =  TMath::Sqrt((dz*dz-d*d)/dq + d*d*h*h/dq/dq)+d*h/dq;
  

  theta = kRaddeg * TMath::ACos(x);
  
  ptrap[0]  = dz/2.;
  ptrap[1]  = theta;
  ptrap[2]  =  0.;
  ptrap[3]  =  d/2;
  ptrap[4]  =  d/x/2;
  ptrap[5]  = ptrap[4];
  ptrap[6]  = 0;
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  ptrap[10] = 0;
  gMC->Gsvolu("BBD1", "TRAP", kSteel, ptrap, 11);
  ptrap[3]  =  d/2-kBBBdd;
  ptrap[4]  = (d/2-kBBBdd)/x;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  gMC->Gsvolu("BBD3", "TRAP", kAir, ptrap, 11);
  gMC->Gspos("BBD3", 1, "BBD1", 0.0, 0.0, 0., 0, "ONLY");
  dx = (kBBMRou + kBBMRin) / 2.;
  dy = ((kBBMRou + kBBMRin)/ 2) *  TMath::Tan(10 * kDegrad) - kBBBdz / 2./ TMath::Cos(10 * kDegrad);
  gMC->Gspos("BBD1", 1, "BBCE", dx, dy,   dz/2. + kBBBdz/2., idrotm[2052], "ONLY");  


  ptrap[0]  = dz/2.;
  ptrap[1]  = -theta;
  ptrap[2]  =  0.;
  ptrap[3]  =  d/2;
  ptrap[4]  =  d/2/x;
  ptrap[5]  = ptrap[4];
  ptrap[6]  = 0;
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  ptrap[10] = 0;
  gMC->Gsvolu("BBD2", "TRAP", kSteel, ptrap, 11);
  ptrap[3]  = d/2-kBBBdd;
  ptrap[4]  = (d/2-kBBBdd)/x;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  gMC->Gsvolu("BBD4", "TRAP", kAir, ptrap, 11);
  gMC->Gspos("BBD4", 1, "BBD2", 0.0, 0.0, 0., 0, "ONLY");
  dx = (kBBMRou + kBBMRin) / 2.;
  dy = ((kBBMRou + kBBMRin)/ 2) *  TMath::Tan(10 * kDegrad) - kBBBdz / 2./ TMath::Cos(10 * kDegrad);
  gMC->Gspos("BBD2", 1, "BBCE", dx, dy,   -dz/2. - kBBBdz/2., idrotm[2052], "ONLY");  


  gMC->Gspos("BBMO", 1, "ALIC", 0., 0., + 376. + kBBMdz / 2. + 0.5, 0, "ONLY");  


}

 

//___________________________________________
void AliFRAMEv2::CreateMaterials()
{
  // Creates the materials
  Float_t epsil, stemax, tmaxfd, deemax, stmin;
  
  epsil  = 1.e-4;     // Tracking precision, 
  stemax = -0.01;     // Maximum displacement for multiple scat 
  tmaxfd = -20.;      // Maximum angle due to field deflection 
  deemax = -.3;       // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();


  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  
  //Air
  
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;

  AliMixture(65, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(5, "AIR$      ", aAir, zAir, dAir,4, wAir);
  AliMaterial(9, "ALU      ", 26.98, 13., 2.7, 8.9, 37.2);

  AliMedium(65, "Stainless Steel", 65, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium( 5, "Air", 5, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium( 9, "Aluminum", 9, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

}

//_____________________________________________________________________________
void AliFRAMEv2::Init()
{
  //
  // Initialise the module after the geometry has been defined
  //
    if(AliLog::GetGlobalDebugLevel()>0) {
	printf("%s: **************************************"
	       " FRAME "
	       "**************************************\n",ClassName());
	printf("\n%s:      Version 2 of FRAME initialised, symmetric FRAME\n\n",ClassName());
	printf("%s: **************************************"
	       " FRAME "
	       "**************************************\n",ClassName());
    }
//
// The reference volume id
    fRefVolumeId = gMC->VolId("BREF");
}

Int_t AliFRAMEv2::IsVersion() const 
{
  // Returns the version of the FRAME (1 if no holes, 0 otherwise) 
    Int_t version = 0;
    if (fHoles == 0) version = 1;
    return version;
}

void AliFRAMEv2::StepManager()
{
//
// Stepmanager of AliFRAMEv2.cxx
// Used for recording of reference tracks entering the spaceframe mother volume
//
  Int_t   copy, id;
  
  //
  // Only charged tracks
  if( !(gMC->TrackCharge()) ) return; 
  //
  // Only tracks entering mother volume
  // 

  id=gMC->CurrentVolID(copy);

  if (id != fRefVolumeId)  return;
  if(!gMC->IsTrackEntering()) return;
  //
  // Add the reference track
  //
  AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
}

  

void AliFRAMEv2::MakeHeatScreen(char* name, Float_t dyP, Int_t rot1, Int_t rot2)
{
    // Heat screen panel
    //
    Int_t *idtmed = fIdtmed->GetArray()-1999;
    const Int_t kAir   = idtmed[2004];
    const Int_t kAlu   = idtmed[2008];

    Float_t dx, dy;
    char mname [128];
    char cname [128];
    char t1name[128];
    char t2name[128];
    char t3name[128];
    char t4name[128];
    char t5name[128];
    
    // 
    Float_t dxP = 98.5 - 3.;
    Float_t dzP =  1.2;
    //
    // Mother volume
    Float_t thshM[3];
    thshM[0]  =  dxP / 2.;
    thshM[1]  =  dyP / 2.;
    thshM[2]  =  dzP / 2.;
    sprintf(mname, "BTSH_%s", name);
    gMC->Gsvolu(mname,  "BOX ", kAir, thshM,  3);
    //
    // Aluminum sheet
    thshM[2] = 0.05;
    sprintf(cname, "BTSHA_%s", name);
    gMC->Gsvolu(cname, "BOX ", kAlu, thshM,  3);
    gMC->Gspos(cname, 1, mname, 0., 0., -0.55, 0);
    //
    // Tubes
    Float_t thshT[3];
    thshT[0] = 0.4;
    thshT[1] = 0.5;
    thshT[2] = (dyP / 2. - 8.);
    //
    sprintf(t1name, "BTSHT1_%s", name);
    gMC->Gsvolu(t1name,  "TUBE", kAlu, thshT,  3);
    dx = - dxP / 2. + 8.;
    gMC->Gspos(t1name, 1, mname,  dx, 0., 0., rot1);
    gMC->Gspos(t1name, 2, mname, -dx, 0., 0., rot1);
    //
    sprintf(t2name, "BTSHT2_%s", name);
    sprintf(t3name, "BTSHT3_%s", name);
    sprintf(t4name, "BTSHT4_%s", name);
    sprintf(t5name, "BTSHT5_%s", name);
    thshT[2] = (thshM[1] - 12.);
    gMC->Gsvolu(t2name,  "TUBE", kAlu, thshT,  3);
    thshT[2] = 7.9/2.;
    gMC->Gsvolu(t3name,  "TUBE", kAlu, thshT,  3);
    thshT[2] = 23.9/2.;
    gMC->Gsvolu(t4name,  "TUBE", kAlu, thshT,  3);
    thshT[2] =  9.0/2.;
    gMC->Gsvolu(t5name,  "TUBE", kAlu, thshT,  3);
    gMC->Gspos(t5name, 0, mname, -dx - 4.,  - (dyP / 2. -  7.5), 0.1, rot2);      

    Int_t sig = 1;
    Int_t ipo = 1;
    for (Int_t i = 0; i < 5; i++) {
	sig *= -1;
	dx += 8.00;
	dy = 4. * sig;
	Float_t dy1 =  - (thshM[1] - 15.5) * sig;
	Float_t dy2 =  - (thshM[1] -  7.5) * sig;
	
	gMC->Gspos(t2name, ipo++, mname, dx, dy, 0.0, rot1);
	dx += 6.9;
	gMC->Gspos(t2name, ipo++, mname, dx, dy, 0.0, rot1);      
	
	gMC->Gspos(t3name, i+1,   mname, dx - 3.45, dy1, 0.0, rot2);      
	gMC->Gspos(t4name, i+1,   mname, dx - 3.45, dy2, 0.0, rot2);      
    }
}



