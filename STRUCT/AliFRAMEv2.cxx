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
Revision 1.3  2001/06/22 12:02:20  morsch
Ring locations matching TRD module positions.

Revision 1.2  2001/05/25 07:59:54  morsch
Initialization print-out in debug mode only.

Revision 1.1  2001/05/11 13:18:05  morsch
C++ version of spaceframe with specs according to Jan Bielski Feb. 2001

*/

////////////////////////////////////////////////
//  space frame class                            /
///////////////////////////////////////////////

#include "AliFRAMEv2.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliMC.h"
#include "TSystem.h"
 
ClassImp(AliFRAMEv2)
 
//_____________________________________________________________________________
AliFRAMEv2::AliFRAMEv2()
{
// Constructor
    SetHoles(0);
}

//_____________________________________________________________________________
AliFRAMEv2::AliFRAMEv2(const char *name, const char *title)
  : AliFRAME(name,title)
{
// Constructor
    SetHoles(0);
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

  Int_t idrotm[2199];
  Int_t *idtmed = fIdtmed->GetArray()-1999;
//
// The Space frame
//
//
  Float_t pbox[3], ptrap[11], ptrd1[4], ptube[3], ppgon[10];
  
  Float_t dx, dy, dz;
  Int_t i, j, jmod;
//
// Constants 
  const Float_t kEps   = 0.01;
  const Int_t kAir   = idtmed[2004];
  const Int_t kSteel = idtmed[2064];
  
  const Float_t krad2deg = 180./TMath::Pi();
  const Float_t kdeg2rad = 1./krad2deg;

  Float_t iFrH   = 114.4;
  Float_t ringH  = 4.;
  Float_t ringW  = 10.;
  Float_t longH  = 5.39;
  Float_t longW  = 6.;  
  Float_t dwl    = 3.14;
  Float_t dwh    = 0.96;

// 
  Float_t dymod[3] = {70., 224., 341.};
//  new ?
//  Float_t dymod[3] = {59.5, 178.5, 341.};

//
// Frame mother volume
//
  ptube[0] = 280.;
//  ptube[1] = 428.2;
  ptube[1] = 430.;
  ptube[2] = 376.;
  gMC->Gsvolu("B077", "TUBE", kAir, ptube, 3);
  gMC->Gspos("B077", 1, "ALIC", 0., 0., 0., 0, "ONLY");
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
  gMC->Gspos("B080", 1, "B077",  281.01, 0., 0., 0, "ONLY");
  gMC->Gspos("B080", 2, "B077", -281.01, 0., 0., 0, "ONLY");

//
// Diagonal bars (1) 
//
  Float_t h, d, dq, x, theta;
  
  h  = (dymod[1]-dymod[0]-2.*dol)*.999;
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
  h  = (2.*dymod[0]-2.*dol)*.999;
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

  AliMatrix(idrotm[2198], -theta+180, 0.0, 90.0, 90.0, 90.-theta, 0.0);
  rtheta = (90.-theta)*kdeg2rad;
  ddx = -ddx0-dol*TMath::Tan(theta2);
  ddz = -dol;
  
  drx = TMath::Cos(rtheta) * ddx +TMath::Sin(rtheta) *ddz+pbox[0];
  drz = -TMath::Sin(rtheta) * ddx +TMath::Cos(rtheta) *ddz-pbox[2];
  gMC->Gspos("B051", 1, "BM49", 
	     drx, 0.0, drz,
	     idrotm[2198], "ONLY");

  AliMatrix(idrotm[2197], -theta, 0.0, 90.0, 90.0, 270.-theta, 0.0);
  rtheta = (270.-theta)*kdeg2rad;
  
  drx =  TMath::Cos(rtheta) * ddx +  TMath::Sin(rtheta) * ddz-pbox[0];
  drz = -TMath::Sin(rtheta) * ddx +  TMath::Cos(rtheta) * ddz+pbox[2];
  gMC->Gspos("B051", 2, "BM49", 
	     drx, 0.0, drz,
	     idrotm[2197], "ONLY");

//
// Diagonal bars (3) 
//
  h  = ((dymod[2]-dymod[1])-2.*dol)*.999;
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
  AliMatrix(idrotm[2121],   0.0, 0.0, 90.0, 130.0, 90.0,  40.0);
  AliMatrix(idrotm[2123], 180.0, 0.0, 90.0, 130.0, 90.0,  40.0);
  AliMatrix(idrotm[2125], 180.0, 0.0, 90.0, 150.0, 90.0, 240.0);
  AliMatrix(idrotm[2130],   0.0, 0.0, 90.0, 150.0, 90.0, 240.0);
  AliMatrix(idrotm[2128],   0.0, 0.0, 90.0, 170.0, 90.0,  80.0);
  AliMatrix(idrotm[2135], 180.0, 0.0, 90.0, 190.0, 90.0, 280.0);
  AliMatrix(idrotm[2141], 180.0, 0.0, 90.0, 170.0, 90.0,  80.0);
  AliMatrix(idrotm[2131],   0.0, 0.0, 90.0, 190.0, 90.0, 280.0);
  AliMatrix(idrotm[2101],   0.0, 0.0, 90.0, 350.0, 90.0, 260.0);
  AliMatrix(idrotm[2106], 180.0, 0.0, 90.0, 350.0, 90.0, 260.0);
  AliMatrix(idrotm[2096], 180.0, 0.0, 90.0,  10.0, 90.0, 100.0);
  AliMatrix(idrotm[2098],   0.0, 0.0, 90.0,  10.0, 90.0, 100.0);
  AliMatrix(idrotm[2108],   0.0, 0.0, 90.0,  30.0, 90.0, 300.0);
  AliMatrix(idrotm[2111], 180.0, 0.0, 90.0,  30.0, 90.0, 300.0);
  AliMatrix(idrotm[2113], 180.0, 0.0, 90.0,  50.0, 90.0, 140.0);
  AliMatrix(idrotm[2118],   0.0, 0.0, 90.0,  50.0, 90.0, 140.0);

  AliMatrix(idrotm[2115], 180.0, 0.0, 90.0, 130.0, 90.0, 220.0);
  AliMatrix(idrotm[2112], 180.0, 0.0, 90.0, 50.0, 90.0, 320.0);
  AliMatrix(idrotm[2124], 180.0, 0.0, 90.0, 150.0, 90.0, 60.0);
  AliMatrix(idrotm[2103], 180.0, 0.0, 90.0, 30.0, 90.0, 120.0);
  AliMatrix(idrotm[2140], 180.0, 0.0, 90.0, 170.0, 90.0, 260.0);
  AliMatrix(idrotm[2133], 180.0, 0.0, 90.0, 190.0, 90.0, 100.0);
  AliMatrix(idrotm[2104], 180.0, 0.0, 90.0, 350.0, 90.0, 80.0);
  AliMatrix(idrotm[2095], 180.0, 0.0, 90.0, 10.0, 90.0, 280.0);
  
  AliMatrix(idrotm[2119], 0.0, 0.0, 90.0, 50.0, 90.0, 320.0);
  AliMatrix(idrotm[2139], 0.0, 0.0, 90.0, 150.0, 90.0, 60.0); 
  AliMatrix(idrotm[2107], 0.0, 0.0, 90.0, 30.0, 90.0, 120.0);
  AliMatrix(idrotm[2099], 0.0, 0.0, 90.0, 10.0, 90.0, 280.0);
  AliMatrix(idrotm[2129], 0.0, 0.0, 90.0, 170.0, 90.0, 260.0);
  AliMatrix(idrotm[2144], 0.0, 0.0, 90.0, 190.0, 90.0, 100.0);
  AliMatrix(idrotm[2100], 0.0, 0.0, 90.0, 350.0, 90.0, 80.0);
  
  Float_t rd =  410.56;
  dz = (dymod[1]+dymod[0])/2.;
  Float_t dz2 =  (dymod[1]+dymod[2])/2.;

//
//  phi = 40
//
  Float_t  phi = 40;
  dx = rd * TMath::Sin(phi*kdeg2rad);
  dy = rd * TMath::Cos(phi*kdeg2rad);
  
  gMC->Gspos("B047", 1, "B076", -dx,  dy,  dz, idrotm[2123], "ONLY");
  gMC->Gspos("B047", 2, "B076", -dx,  dy, -dz, idrotm[2121], "ONLY");
  gMC->Gspos("B047", 3, "B076",  dx,  dy,  dz, idrotm[2113], "ONLY");
  gMC->Gspos("B047", 4, "B076",  dx,  dy, -dz, idrotm[2118], "ONLY");

  gMC->Gspos("B045", 1, "B076", -dx,  dy,  dz2, idrotm[2115], "ONLY");
  gMC->Gspos("B045", 2, "B076", -dx,  dy, -dz2, idrotm[2121], "ONLY"); // ?
  gMC->Gspos("B045", 3, "B076",  dx,  dy,  dz2, idrotm[2112], "ONLY");
  gMC->Gspos("B045", 4, "B076",  dx,  dy, -dz2, idrotm[2119], "ONLY");

  gMC->Gspos("BM49", 1, "B076",  dx,  dy,  0., idrotm[2112], "ONLY");
  gMC->Gspos("BM49", 2, "B076", -dx,  dy,  0., idrotm[2115], "ONLY");

//
//  phi = 60
//

  phi = 60;
  dx = rd * TMath::Sin(phi*kdeg2rad);
  dy = rd * TMath::Cos(phi*kdeg2rad);
  gMC->Gspos("B047", 5, "B076", -dx,  dy,  dz, idrotm[2125], "ONLY");
  gMC->Gspos("B047", 6, "B076", -dx,  dy, -dz, idrotm[2130], "ONLY");
  gMC->Gspos("B047", 7, "B076",  dx,  dy,  dz, idrotm[2111], "ONLY");
  gMC->Gspos("B047", 8, "B076",  dx,  dy, -dz, idrotm[2108], "ONLY");

  gMC->Gspos("B045", 5, "B076", -dx,  dy,  dz2, idrotm[2124], "ONLY");
  gMC->Gspos("B045", 6, "B076", -dx,  dy, -dz2, idrotm[2139], "ONLY");
  gMC->Gspos("B045", 7, "B076",  dx,  dy,  dz2, idrotm[2103], "ONLY");
  gMC->Gspos("B045", 8, "B076",  dx,  dy, -dz2, idrotm[2107], "ONLY");

  gMC->Gspos("BM49", 3, "B076",  dx,  dy,  0., idrotm[2103], "ONLY");
  gMC->Gspos("BM49", 4, "B076", -dx,  dy,  0., idrotm[2124], "ONLY");
//
//  phi = 80
//

  phi = 80;
  dx = rd * TMath::Sin(phi*kdeg2rad);
  dy = rd * TMath::Cos(phi*kdeg2rad);
  gMC->Gspos("B047",  9, "B076", -dx,  dy,  dz, idrotm[2141], "ONLY");
  gMC->Gspos("B047", 10, "B076", -dx,  dy, -dz, idrotm[2128], "ONLY");
  gMC->Gspos("B047", 11, "B076",  dx,  dy,  dz, idrotm[2096], "ONLY");
  gMC->Gspos("B047", 12, "B076",  dx,  dy, -dz, idrotm[2098], "ONLY");

  gMC->Gspos("B047", 13, "B076", -dx, -dy,  dz, idrotm[2135], "ONLY");
  gMC->Gspos("B047", 14, "B076", -dx, -dy, -dz, idrotm[2131], "ONLY");
  gMC->Gspos("B047", 15, "B076",  dx, -dy,  dz, idrotm[2106], "ONLY");
  gMC->Gspos("B047", 16, "B076",  dx, -dy, -dz, idrotm[2101], "ONLY");

  gMC->Gspos("B045",  9, "B076", -dx,  dy,  dz2, idrotm[2140], "ONLY");
  gMC->Gspos("B045", 10, "B076", -dx,  dy, -dz2, idrotm[2129], "ONLY");
  gMC->Gspos("B045", 11, "B076",  dx,  dy,  dz2, idrotm[2095], "ONLY");
  gMC->Gspos("B045", 12, "B076",  dx,  dy, -dz2, idrotm[2099], "ONLY");

  gMC->Gspos("B045", 13, "B076", -dx, -dy,  dz2, idrotm[2133], "ONLY");
  gMC->Gspos("B045", 14, "B076", -dx, -dy, -dz2, idrotm[2144], "ONLY");
  gMC->Gspos("B045", 15, "B076",  dx, -dy,  dz2, idrotm[2104], "ONLY");
  gMC->Gspos("B045", 16, "B076",  dx, -dy, -dz2, idrotm[2100], "ONLY");

  gMC->Gspos("BM49", 5, "B076",  dx,  dy,  0., idrotm[2095], "ONLY");
  gMC->Gspos("BM49", 6, "B076", -dx,  dy,  0., idrotm[2140], "ONLY");
  gMC->Gspos("BM49", 7, "B076",  dx, -dy,  0., idrotm[2104], "ONLY");
  gMC->Gspos("BM49", 8, "B076", -dx, -dy,  0., idrotm[2133], "ONLY");


// The internal frame
//
  char*  module[3] = {"B071\0", "B074\0", "B075\0"};
//
//
//  Mother Volumes
//
  ptrd1[0] = 50.18;
  ptrd1[1] = 70.35;
//  ptrd1[0] = 50.10;
//  ptrd1[1] = 70.20;

  ptrd1[2] = 375.5;
  ptrd1[3] =  57.2;  
  for (jmod = 0; jmod < 3; jmod++)
  {
      gMC->Gsvolu(module[jmod], "TRD1", kAir, ptrd1, 4);
  }

  Int_t mod[18] = {1, 1, 1, 0, 0, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 1, 1};
  Int_t rot1[18] = {53, 52, 49, 69, 30, 71, 64, 83, 81, 151, 150, 145, 146, 147, 148, 149, 60, 51};
  Int_t rot2[18] = {23, 22, 21, 24, 14, 20, 136, 137, 36, 28, 67, 25, 50, 26, 27, 29, 138, 2};
  
  
  Float_t r      = 341.8;
  Float_t rout1  = 410.564;
  Float_t rout2  = 415.2;
  Int_t modcount[3] = {0, 0, 0};
  
  for (i=0; i<18; i++) {
      Float_t phi = i*20.;
      Float_t phi2 = 270+phi;
      if (phi2 >= 360.) phi2-=360.;
      
      dx =  TMath::Sin(phi*kdeg2rad)*r;
      dy = -TMath::Cos(phi*kdeg2rad)*r;
      modcount[mod[i]]++;
      
      AliMatrix(idrotm[2000+rot1[i]],  90.0, phi, 0., 0., 90., phi2);  
      gMC->Gspos(module[mod[i]], modcount[mod[i]], "B077", dx, dy, 0., idrotm[2000+rot1[i]], "ONLY");
//
//    Position elements of outer Frame
//
      dx =  TMath::Sin(phi*kdeg2rad)*rout1;
      dy = -TMath::Cos(phi*kdeg2rad)*rout1;
      for (j = 0; j < 3; j++)
      {
	  dz = dymod[j];
	  gMC->Gspos("B042", 6*i+2*j+1, "B076", dx, dy,  dz, idrotm[2000+rot1[i]], "ONLY");	  
	  gMC->Gspos("B042", 6*i+2*j+2, "B076", dx, dy, -dz, idrotm[2000+rot1[i]], "ONLY");	  
      }

      phi = i*20.+10;
      phi2 = 270+phi;
      AliMatrix(idrotm[2000+rot2[i]],  90.0, phi, 90., phi2, 0., 0.);  

      dx =  TMath::Sin(phi*kdeg2rad)*rout2;
      dy = -TMath::Cos(phi*kdeg2rad)*rout2;
      gMC->Gspos("B033", i+1, "B076", dx, dy,  0., idrotm[2000+rot2[i]], "ONLY");	  
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
  
  pbox[0] = dwh;
  pbox[1] = dwl;
  pbox[2] = (iFrH-ringH-longH)/2.;
  gMC->Gsvolu("B063", "BOX ", kSteel, pbox, 3);
  pbox[0] = dwh-0.2;
  pbox[1] = dwl-0.4;
  
  gMC->Gsvolu("B064", "BOX ", kAir, pbox, 3);
  gMC->Gspos("B064", 1, "B063", 0.2, 0., 0., 0, "ONLY");
  
  AliMatrix(idrotm[2013],  90.0,   0.0,  90.0, 270.0,   0.0,   0.0);  
  AliMatrix(idrotm[2001], 100.0,   0.0,  90.0, 270.0,  10.0,   0.0);
  AliMatrix(idrotm[2003], 100.0,   0.0,  90.0,  90.0,  10.0,   0.0);
  AliMatrix(idrotm[2011], 100.0, 180.0,  90.0, 270.0,  10.0, 180.0);
  AliMatrix(idrotm[2017], 100.0, 180.0,  90.0,  90.0,  10.0, 180.0);
  AliMatrix(idrotm[2006],  10.0,   0.0,  80.0, 180.0,  90.0,  90.0);
  AliMatrix(idrotm[2010],  10.0,   0.0,  80.0, 180.0,  90.0, 270.0);
  AliMatrix(idrotm[2012],  10.0, 180.0,  80.0,   0.0,  90.0,  90.0);
  AliMatrix(idrotm[2018],  10.0, 180.0,  80.0,   0.0,  90.0, 270.0);
  AliMatrix(idrotm[2007], 170.0, 180.0,  80.0, 180.0,  90.0,  90.0);
  AliMatrix(idrotm[2009], 170.0, 180.0,  80.0, 180.0,  90.0, 270.0);
  AliMatrix(idrotm[2016], 170.0,   0.0,  80.0,   0.0,  90.0,  90.0);
  AliMatrix(idrotm[2015], 170.0,   0.0,  80.0,   0.0,  90.0, 270.0);
  AliMatrix(idrotm[2004], 170.0,   0.0,  90.0,  90.0,  80.0,   0.0);
  AliMatrix(idrotm[2005], 170.0, 180.0,  90.0,  90.0,  80.0, 180.0);
  AliMatrix(idrotm[2019],  90.0, 180.0,  90.0,  90.0,   0.0,   0.0);
//
// web frame diagonal (outer)
//  
  h  = 106.2;
  d  = 2.*dwl;
  dz = dymod[2]-dymod[1]-2.*dwl;
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
  dz = dymod[1]-dymod[0]-2.*dwl;
  dq = h*h+dz*dz;
   x  =  TMath::Sqrt((dz*dz-d*d)/dq + d*d*h*h/dq/dq)+d*h/dq;
  

  theta = krad2deg * TMath::ACos(x);
  
  ptrap[0]  = (dymod[1]-dymod[0]-2.*dwl)/2.;
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
  
  for (i = 0; i< 3; i++)
  {
// ring bars
      for (jmod = 0; jmod<3; jmod++) {
	  gMC->Gspos("B072", 6*i+jmod+1, module[jmod], 0,  dymod[i], dz, 0, "ONLY");
	  gMC->Gspos("B072", 6*i+jmod+4, module[jmod], 0, -dymod[i], dz, idrotm[2013], "ONLY");      

// 0-deg web
	  gMC->Gspos("B063", 12*i+jmod+1,  module[jmod],  60.0732,  dymod[i], 4.6669, idrotm[2003], "ONLY");
	  gMC->Gspos("B063", 12*i+jmod+4,  module[jmod],  60.0732, -dymod[i], 4.6669, idrotm[2001], "ONLY");      
	  gMC->Gspos("B063", 12*i+jmod+7,  module[jmod], -60.0732,  dymod[i], 4.6669, idrotm[2017], "ONLY");
	  gMC->Gspos("B063", 12*i+jmod+10, module[jmod], -60.0732, -dymod[i], 4.6669, idrotm[2011], "ONLY");      
      }
  }
  
// outer diagonal web

  dy = (dymod[2]+dymod[1])/2.;
  for (jmod = 0; jmod<3; jmod++) {
      gMC->Gspos("B065", 4*jmod+1, module[jmod],  60.0732,   dy, 4.6669, idrotm[2010], "ONLY");
      gMC->Gspos("B065", 4*jmod+2, module[jmod],  60.0732,  -dy, 4.6669, idrotm[2006], "ONLY");
      gMC->Gspos("B065", 4*jmod+3, module[jmod], -60.0732,   dy, 4.6669, idrotm[2018], "ONLY");
      gMC->Gspos("B065", 4*jmod+4, module[jmod], -60.0732,  -dy, 4.6669, idrotm[2012], "ONLY");
  }
  

  dy = (dymod[1]+dymod[0])/2.;

  for (jmod = 0; jmod<3; jmod++) {
      gMC->Gspos("B067", 4*jmod+1, module[jmod],  60.0732,   dy, 4.6669, idrotm[2009], "ONLY");
      gMC->Gspos("B067", 4*jmod+2, module[jmod],  60.0732,  -dy, 4.6669, idrotm[2007], "ONLY");
      gMC->Gspos("B067", 4*jmod+3, module[jmod], -60.0732,   dy, 4.6669, idrotm[2015], "ONLY");
      gMC->Gspos("B067", 4*jmod+4, module[jmod], -60.0732,  -dy, 4.6669, idrotm[2016], "ONLY");
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
  gMC->Gsvolu("B059", "TRAP", kSteel, ptrap, 11);
  ptrap[0]  =  2.2;
  ptrap[4]  =  2.15;
  ptrap[5]  = ptrap[4];
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  gMC->Gsvolu("B062", "TRAP", kAir, ptrap, 11);
  gMC->Gspos("B062", 1, "B059", 0.0, -0.3, 0., 0, "ONLY");
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

  for (jmod = 0; jmod<3; jmod++) {
      gMC->Gspos("BA59", 2*jmod+1, module[jmod],  49.6476, 0.0, dz, idrotm[2005], "ONLY");
      gMC->Gspos("BA59", 2*jmod+2, module[jmod], -49.6476, 0.0, dz, idrotm[2004], "ONLY");
  }
//
//
// Spacer between TRD and TOF modules
  ptrap[0]  =   0.99;
  ptrap[1]  =  10.00;
  ptrap[2]  =   0.00;
  ptrap[3]  = 350.00;
  ptrap[4]  =    0.9;
  ptrap[5]  = ptrap[4];
  ptrap[6]  = 0;
  ptrap[7]  = ptrap[3];
  ptrap[8]  = ptrap[4];
  ptrap[9]  = ptrap[4];
  ptrap[10] = 0;
  gMC->Gsvolu("B056", "TRAP", kSteel, ptrap, 11);
  for (jmod = 0; jmod<3; jmod++) {
      gMC->Gspos("B056", 2*jmod+1, module[jmod],  61.9607, 0.0, 27.2, 0           , "ONLY");
      gMC->Gspos("B056", 2*jmod+2, module[jmod], -61.9607, 0.0, 27.2, idrotm[2019], "ONLY");
  }
//
// Mother volumes for TRD and TOF 
// 
  if (!fHoles) {
      
      ptrd1[0] = 49.8065;
      ptrd1[1] = 62.8535;
      ptrd1[2] = 375.5;
      ptrd1[3] = 37.;
      gMC->Gsvolu("BTR1", "TRD1", kAir, ptrd1, 4);
      gMC->Gsvolu("BTR2", "TRD1", kAir, ptrd1, 4);
      gMC->Gsvolu("BTR3", "TRD1", kAir, ptrd1, 4);  
      
      ptrd1[0] = 63.2061;
      ptrd1[1] = 68.3192;
      ptrd1[2] = 375.5;
      ptrd1[3] = 14.5;
      gMC->Gsvolu("BTO1", "TRD1", kAir, ptrd1, 4);
      gMC->Gsvolu("BTO2", "TRD1", kAir, ptrd1, 4);
      gMC->Gsvolu("BTO3", "TRD1", kAir, ptrd1, 4);  
      
      gMC->Gspos("BTR1", 1, "B071", 0., 0., -10.8, 0, "ONLY");
      gMC->Gspos("BTR2", 1, "B074", 0., 0., -10.8, 0, "ONLY");
      gMC->Gspos("BTR3", 1, "B075", 0., 0., -10.8, 0, "ONLY");
      
      gMC->Gspos("BTO1", 1, "B071", 0., 0.,  42.69, 0, "ONLY");
      gMC->Gspos("BTO2", 1, "B074", 0., 0.,  42.69, 0, "ONLY");
      gMC->Gspos("BTO3", 1, "B075", 0., 0.,  42.69, 0, "ONLY");
  } else {
      ptrd1[0] = 49.8065;
      ptrd1[1] = 62.8535;
      ptrd1[2] = 375.5;
      ptrd1[3] = 37;
      gMC->Gsvolu("BTR1", "TRD1", kAir, ptrd1, 4);
      ptrd1[2] = 156.75;
      gMC->Gsvolu("BTR2", "TRD1", kAir, ptrd1, 4);
      ptrd1[2] =  79.75;
      gMC->Gsvolu("BTR3", "TRD1", kAir, ptrd1, 4);  

      ptrd1[0] = 63.2061;
      ptrd1[1] = 68.3192;
      ptrd1[2] = 375.5;
      ptrd1[3] = 14.5;
      gMC->Gsvolu("BTO1", "TRD1", kAir, ptrd1, 4);
      ptrd1[2] = 156.75;
      gMC->Gsvolu("BTO2", "TRD1", kAir, ptrd1, 4);
      ptrd1[2] =  79.75;
      gMC->Gsvolu("BTO3", "TRD1", kAir, ptrd1, 4);  
      
      gMC->Gspos("BTR1", 1, "B071", 0., 0., -10.8, 0, "ONLY");

      gMC->Gspos("BTR2", 1, "B074", 0., -218.75, -10.8, idrotm[152], "ONLY");
      gMC->Gspos("BTR2", 2, "B074", 0.,  218.75, -10.8,           0, "ONLY");

      gMC->Gspos("BTR3", 1, "B075", 0., -295.75, -10.8, idrotm[152], "ONLY");
      gMC->Gspos("BTR3", 2, "B075", 0.,  295.75, -10.8,           0, "ONLY");


      gMC->Gspos("BTO1", 1, "B071", 0., 0., -42.7, 0, "ONLY");

      gMC->Gspos("BTO2", 1, "B074", 0., -218.75, -42.7, idrotm[152], "ONLY");
      gMC->Gspos("BTO2", 2, "B074", 0.,  218.75, -42.7,           0, "ONLY");

      gMC->Gspos("BTO3", 1, "B075", 0., -295.75, -42.7, idrotm[152], "ONLY");
      gMC->Gspos("BTO3", 2, "B075", 0.,  295.75, -42.7,           0, "ONLY");
      
  }
}
 

//___________________________________________
void AliFRAMEv2::CreateMaterials()
{

  Float_t epsil, stemax, tmaxfd, deemax, stmin;
  
  epsil  = 1.e-4;     // Tracking precision, 
  stemax = -0.01;   // Maximum displacement for multiple scat 
  tmaxfd = -20.;    // Maximum angle due to field deflection 
  deemax = -.3;     // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };

  AliMixture(65, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  AliMaterial(5, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500.);
  AliMedium(65, "Stainless Steel", 65, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium( 5, "Air            ", 5, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
}

//_____________________________________________________________________________
void AliFRAMEv2::Init()
{
  //
  // Initialise the module after the geometry has been defined
  //
    if(fDebug) {
	printf("%s: **************************************"
	       " FRAME "
	       "**************************************\n",ClassName());
	printf("\n%s:      Version 2 of FRAME initialised, symmetric FRAME\n\n",ClassName());
	printf("%s: **************************************"
	       " FRAME "
	       "**************************************\n",ClassName());
    }
}








