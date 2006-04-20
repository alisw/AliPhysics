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
Revision 1.4  2006/04/16 22:29:05  hristov
Coding conventions (Annalisa)

Revision 1.3  2006/03/12 14:38:05  arcelli
 Changes for TOF Reconstruction using TGeo

Revision 1.2  2006/02/28 10:38:00  decaro
AliTOFGeometry::fAngles, AliTOFGeometry::fHeights, AliTOFGeometry::fDistances arrays: dimension definition in the right location

Revision 1.1  2005/12/15 08:55:33  decaro
New TOF geometry description (V5) -G. Cara Romeo and A. De Caro

Revision 0.1  2005/07/19 G. Cara Romeo and A. De Caro
        Modify Global methods IsInsideThePad & DistanceToPad
               according to the new TOF geometry
        Implement Global  methods GetPadDx & GetPadDy & GetPadDz
        Implement Private methods Translation & Rotation & InverseRotation
        Modify Global methods GetDetID & GetPlate & GetSector &
                              GetStrip & GetPadX & GetPadZ
               according to the new TOF geometry
        Modify Global methods GetPos & GetX & GetY & GetZ
               according to the new TOF geometry
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TOF Geometry class (new version)                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TGeoManager.h"

#include "AliConst.h"
#include "AliLog.h"

#include "AliTOFGeometryV5.h"

extern TGeoManager *gGeoManager;

ClassImp(AliTOFGeometryV5)


const Float_t AliTOFGeometryV5::fgkZlenA    = 370.6*2.; // length (cm) of the A module
const Float_t AliTOFGeometryV5::fgkZlenB    = 146.5;    // length (cm) of the B module
const Float_t AliTOFGeometryV5::fgkZlenC    = 170.45;   // length (cm) of the C module
const Float_t AliTOFGeometryV5::fgkMaxhZtof = 370.6;    // Max half z-size of TOF (cm)

const Float_t AliTOFGeometryV5::fgkxTOF     = 371.-0.01;// Inner radius of the TOF for Reconstruction (cm)
const Float_t AliTOFGeometryV5::fgkRmin     = 370.-0.01;// Inner radius of the TOF (cm)
const Float_t AliTOFGeometryV5::fgkRmax     = 399.-0.01;// Outer radius of the TOF (cm)

//_____________________________________________________________________________
AliTOFGeometryV5::AliTOFGeometryV5()
  :AliTOFGeometry()
{
  //
  // AliTOFGeometryV5 default constructor
  //

  AliTOFGeometry::fNStripC     = kNStripC;       // number of strips in C type module

  AliTOFGeometry::fZlenA       = fgkZlenA;       // length of the TOF supermodule (cm)
  AliTOFGeometry::fZlenB       = fgkZlenB;       // length of the B module (cm)
  AliTOFGeometry::fZlenC       = fgkZlenC;       // length of the C module (cm)
  AliTOFGeometry::fMaxhZtof    = fgkMaxhZtof;    // Max half z-size of TOF supermodule (cm)

  AliTOFGeometry::fxTOF   = fgkxTOF;           // Inner radius of the TOF for Reconstruction (cm)
  AliTOFGeometry::fRmin   = fgkRmin;           // Inner radius of the TOF (cm)
  AliTOFGeometry::fRmax   = fgkRmax;           // Outer radius of the TOF (cm)

  Init();

}

//_____________________________________________________________________________
AliTOFGeometryV5::~AliTOFGeometryV5()
{
  //
  // AliTOFGeometryV5 destructor
  //

}
//_____________________________________________________________________________
void AliTOFGeometryV5::ImportGeometry(){
  TGeoManager::Import("geometry.root");
}
//_____________________________________________________________________________
void AliTOFGeometryV5::Init()
{
  //
  // Initialize strip Tilt Angles, Heights and Distances
  //
  // Strips Tilt Angles
 
  // For each strip to be positoned in FLTA/FLTB/FLTC,
  // define 3 arrays containing:
  //   the angle of the normal with respect to the Y axis of FLTA/FLTB/FLTC
  //   the Y of the center with respect to the FLTA/FLTB/FLTC reference frame
  //   the Z of the center with respect to the BT01/BT02/BT03 reference frame


  fPhiSec   = 360./kNSectors;

  Float_t const kangles[kNPlates][kMaxNstrip] ={
    { 43.99,  43.20,  42.40,  41.59,  40.77,  39.94,  39.11,  38.25,  37.40,  36.53,
      35.65,  34.76,  33.87,  32.96,  32.05,  31.13,  30.19,  29.24,  12.33,  0.00},

    { 27.26,  26.28,  25.30,  24.31,  23.31,  22.31,  21.30,  20.29,  19.26,  18.24,
      17.20,  16.16,  15.11,  14.05,  13.00,  11.93,  10.87,   9.80,   8.74,  0.00},

    {  0.00,   6.30,   5.31,   4.25,   3.19,   2.12,   1.06,   0.00,  -1.06,  -2.12,
      -3.19,  -4.25,  -5.31,  -6.30,   0.00,   0.00,   0.00,   0.00,   0.00,  0.00},

    { -8.74,  -9.80, -10.87, -11.93, -13.00, -14.05, -15.11, -16.16, -17.20, -18.24,
     -19.26, -20.29, -21.30, -22.31, -23.31, -24.31, -25.30, -26.28, -27.26,  0.00},
    
    {-12.33, -29.24, -30.19, -31.13, -32.05, -32.96, -33.87, -34.76, -35.65, -36.53,
     -37.40, -38.25, -39.11, -39.94, -40.77, -41.59, -42.40, -43.20, -43.99,  0.00}
  };

  Float_t const kheights[kNPlates][kMaxNstrip]= {
    {-8.2,  -7.5,  -8.2,  -7.7,  -8.1,  -7.6,  -7.7,  -7.7,  -7.7,  -7.7,
     -7.5,  -7.2,  -7.3,  -7.5,  -7.6,  -7.8,  -8.3,  -9.3,  -3.1,   0.0},

    {-7.9,  -8.1,  -8.5,  -9.0, -10.1,  -3.9,  -5.9,  -7.7, -10.1,  -3.6,
     -5.8,  -8.0, -10.4,  -4.4,  -7.2, -10.2,  -4.6,  -7.4, -10.4,   0.0},

    {-2.5, -10.4,  -5.0,  -9.9,  -4.8,  -9.9,  -4.7, -10.2,  -4.7,  -9.9,
     -4.8,  -9.9,  -5.0, -10.4,  -2.5,   0.0,   0.0,   0.0,   0.0,   0.0},

    {-10.4, -7.4,  -4.6, -10.2,  -7.2,  -4.4, -10.4,  -8.0,  -5.8,  -3.6,
     -10.1,  -7.7, -5.9,  -3.9, -10.1,  -9.0,  -8.5,  -8.1,  -7.9,   0.0},

    { -3.1,  -9.3, -8.3,  -7.8,  -7.6,  -7.5,  -7.3,  -7.2,  -7.5,  -7.7,
      -7.7,  -7.7, -7.7,  -7.6,  -8.1,  -7.7,  -8.2,  -7.5,  -8.2,   0.0}
  };


  Float_t const kdistances[kNPlates][kMaxNstrip]= {
    { 364.1,  354.9,  344.5,  335.4,  325.5,  316.6,  307.2,  298.0,  288.9,  280.0,
      271.3,  262.7,  254.0,  244.8,  236.1,  227.7,  219.1,  210.3,  205.7,    0.0},

    { 194.2,  186.1,  177.9,  169.8,  161.5,  156.3,  147.8,  139.4,  130.9,  125.6,
      117.3,  109.2,  101.1,   95.3,   87.1,   79.2,   73.0,   65.1,   57.6,    0.0},

    {  49.5,   41.3,   35.3,   27.8,   21.2,   13.9,    7.0,    0.0,   -7.0,  -13.9,
      -21.2,  -27.8,  -35.3,  -41.3,  -49.5,    0.0,    0.0,    0.0,    0.0,    0.0},

    { -57.6,  -65.1,  -73.0,  -79.2,  -87.1,  -95.3, -101.1, -109.2, -117.3, -125.6,
     -130.9, -139.4, -147.8, -156.3, -161.5, -169.8, -177.9, -186.1, -194.2,    0.0},

    {-205.7, -210.3, -219.1, -227.7, -236.1, -244.8, -254.0, -262.7, -271.3, -280.0,
     -288.9, -298.0, -307.2, -316.6, -325.5, -335.4, -344.5, -354.9, -364.1,    0.0}
  };


  for (Int_t iplate = 0; iplate < kNPlates; iplate++) {
    for (Int_t istrip = 0; istrip < kMaxNstrip; istrip++) {
      AliTOFGeometry::fAngles[iplate][istrip]   = kangles[iplate][istrip];
      AliTOFGeometry::fHeights[iplate][istrip]  = kheights[iplate][istrip];
      AliTOFGeometry::fDistances[iplate][istrip]= kdistances[iplate][istrip];
    }
  }

}

//_____________________________________________________________________________
Float_t AliTOFGeometryV5::DistanceToPadPar(Int_t *det, Float_t *pos, Float_t *dist3d) const
{
//
// Returns distance of  space point with coor pos (x,y,z) (cm) wrt 
// pad with Detector Indices idet (iSect,iPlate,iStrip,iPadX,iPadZ) 
//
    
  //Transform pos into Sector Frame

  Float_t x = pos[0];
  Float_t y = pos[1];
  Float_t z = pos[2];

  Float_t radius = TMath::Sqrt(x*x+y*y);
  //Float_t phi=TMath::ATan(y/x);	
  //if(phi<0) phi = k2PI+phi; //2.*TMath::Pi()+phi;
  Float_t phi = TMath::Pi()+TMath::ATan2(-y,-x);	
  //  Get the local angle in the sector philoc
  Float_t angle   = phi*kRaddeg-( Int_t (kRaddeg*phi/fPhiSec) + 0.5)*fPhiSec;
  Float_t xs = radius*TMath::Cos(angle/kRaddeg);
  Float_t ys = radius*TMath::Sin(angle/kRaddeg);
  Float_t zs = z;

  // Do the same for the selected pad

  Float_t g[3];
  GetPosPar(det,g);

  Float_t padRadius = TMath::Sqrt(g[0]*g[0]+g[1]*g[1]);
  //Float_t padPhi = TMath::ATan(g[1]/g[0]);	
  //if(padPhi<0) padPhi = k2Pi + padPhi;
  Float_t padPhi = TMath::Pi()+TMath::ATan2(-g[1],-g[0]);	

  //  Get the local angle in the sector philoc
  Float_t padAngle = padPhi*kRaddeg-( Int_t (padPhi*kRaddeg/fPhiSec)+ 0.5) * fPhiSec;
  Float_t padxs = padRadius*TMath::Cos(padAngle/kRaddeg);
  Float_t padys = padRadius*TMath::Sin(padAngle/kRaddeg);
  Float_t padzs = g[2];
  
  //Now move to local pad coordinate frame. Translate:
  
  Float_t xt = xs-padxs;
  Float_t yt = ys-padys;
  Float_t zt = zs-padzs;
  //Now Rotate:
  
  Float_t alpha = GetAngles(det[1],det[2]);
  Float_t xr =  xt*TMath::Cos(alpha/kRaddeg)+zt*TMath::Sin(alpha/kRaddeg);
  Float_t yr =  yt;
  Float_t zr = -xt*TMath::Sin(alpha/kRaddeg)+zt*TMath::Cos(alpha/kRaddeg);

  Float_t dist = TMath::Sqrt(xr*xr+yr*yr+zr*zr);

  if (dist3d){
    dist3d[0] = xr;
    dist3d[1] = yr;
    dist3d[2] = zr;
  }

  return dist;

}

//_____________________________________________________________________________
Bool_t AliTOFGeometryV5::IsInsideThePadPar(Int_t *det, Float_t *pos) const
{
//
// Returns true if space point with coor pos (x,y,z) (cm) falls 
// inside pad with Detector Indices idet (iSect,iPlate,iStrip,iPadX,iPadZ) 
//

  Bool_t isInside=false; 

  /*
  const Float_t khhony    = 1.0          ; // heigth of HONY  Layer
  const Float_t khpcby    = 0.08         ; // heigth of PCB   Layer
  const Float_t khrgly    = 0.055        ; // heigth of RED GLASS  Layer
  const Float_t khglfy    = 0.285        ; // heigth of GLASS+FISHLINE  Layer
  const Float_t khcpcby   = 0.16         ; // heigth of PCB  Central Layer
  //const Float_t kwcpcbz   = 12.4         ; // z dimension of PCB  Central Layer
  const Float_t khstripy = 2.*khhony+2.*khpcby+4.*khrgly+2.*khglfy+khcpcby;//3.11
  //const Float_t kwstripz = kwcpcbz;
  //const Float_t klstripx = fgkStripLength;
  */

  const Float_t khsensmy = 0.5;//0.05;//0.11;//0.16;//          // heigth of Sensitive Layer

  //Transform pos into Sector Frame

  Float_t x = pos[0];
  Float_t y = pos[1];
  Float_t z = pos[2];

  Float_t radius = TMath::Sqrt(x*x+y*y);
  Float_t phi = TMath::Pi()+TMath::ATan2(-y,-x);	

  //  Get the local angle in the sector philoc
  Float_t angle = phi*kRaddeg-( Int_t (kRaddeg*phi/fPhiSec) + 0.5) *fPhiSec;
  Float_t xs = radius*TMath::Cos(angle/kRaddeg);
  Float_t ys = radius*TMath::Sin(angle/kRaddeg);
  Float_t zs = z;

  // Do the same for the selected pad

  Float_t g[3];
  GetPosPar(det,g);

  Float_t padRadius = TMath::Sqrt(g[0]*g[0]+g[1]*g[1]);
  Float_t padPhi = TMath::Pi()+TMath::ATan2(-g[1],-g[0]);	

  //  Get the local angle in the sector philoc
  Float_t padAngle = padPhi*kRaddeg-( Int_t (padPhi*kRaddeg/fPhiSec)+ 0.5) * fPhiSec; 
  Float_t padxs = padRadius*TMath::Cos(padAngle/kRaddeg);
  Float_t padys = padRadius*TMath::Sin(padAngle/kRaddeg);
  Float_t padzs = g[2];

  //Now move to local pad coordinate frame. Translate:

  Float_t xt = xs-padxs;
  Float_t yt = ys-padys;
  Float_t zt = zs-padzs;

  //Now Rotate:

  Float_t alpha = GetAngles(det[1],det[2]);
  Float_t xr =  xt*TMath::Cos(alpha/kRaddeg)+zt*TMath::Sin(alpha/kRaddeg);
  Float_t yr =  yt;
  Float_t zr = -xt*TMath::Sin(alpha/kRaddeg)+zt*TMath::Cos(alpha/kRaddeg);

  if(TMath::Abs(xr)<=khsensmy*0.5 && TMath::Abs(yr)<= (fgkXPad*0.5) && TMath::Abs(zr)<= (fgkZPad*0.5))
    isInside=true;
  return isInside;

}


//_____________________________________________________________________________
Float_t AliTOFGeometryV5::DistanceToPad(Int_t *det, TGeoHMatrix mat, Float_t *pos, Float_t *dist3d) const
{
//
// Returns distance of  space point with coor pos (x,y,z) (cm) wrt 
// pad with Detector Indices idet (iSect,iPlate,iStrip,iPadX,iPadZ) 
//
  if (!gGeoManager) {
    printf("ERROR: no TGeo\n");
    return 0.;
  }
  Double_t vecg[3];
  vecg[0]=pos[0];
  vecg[1]=pos[1];
  vecg[2]=pos[2];
  Double_t veclr[3]={-1.,-1.,-1.};
  Double_t vecl[3]={-1.,-1.,-1.};
  mat.MasterToLocal(vecg,veclr);  
  vecl[0]=veclr[1];
  vecl[1]=veclr[0];
  //take into account reflections 
  if(det[1]>-1)vecl[2]=-veclr[2];

  Float_t dist = TMath::Sqrt(vecl[0]*vecl[0]+vecl[1]*vecl[1]+vecl[2]*vecl[2]);


  if (dist3d){
    dist3d[0] = vecl[0];
    dist3d[1] = vecl[1];
    dist3d[2] = vecl[2];
  }

  return dist;

}


//_____________________________________________________________________________
Bool_t AliTOFGeometryV5::IsInsideThePad( Int_t *det, TGeoHMatrix mat, Float_t *pos) const
{
//
// Returns true if space point with coor pos (x,y,z) (cm) falls 
// inside pad with Detector Indices idet (iSect,iPlate,iStrip,iPadX,iPadZ) 
//

  const Float_t khsensmy = 0.5;      // heigth of Sensitive Layer
  Double_t vecg[3];
  vecg[0]=pos[0];
  vecg[1]=pos[1];
  vecg[2]=pos[2];
  Double_t veclr[3]={-1.,-1.,-1.};
  Double_t vecl[3]={-1.,-1.,-1.};
  mat.MasterToLocal(vecg,vecl);  
  vecl[0]=veclr[1];
  vecl[1]=veclr[0];
  //take into account reflections 
  if(det[1]>-1)vecl[2]=-veclr[2];

  Float_t xr = vecl[0];
  Float_t yr = vecl[1];
  Float_t zr = vecl[2];

  Bool_t isInside=false; 
  if(TMath::Abs(xr)<= khsensmy*0.5 && TMath::Abs(yr)<= (fgkXPad*0.5) && TMath::Abs(zr)<= (fgkZPad*0.5))
    isInside=true; 
  return isInside;

}
//_____________________________________________________________________________
//_____________________________________________________________________________
Float_t AliTOFGeometryV5::GetX(Int_t *det) const
{
  //
  // Returns X coordinate (cm)
  //

  Int_t isector = det[0];
  Int_t iplate  = det[1];
  Int_t istrip  = det[2];
  Int_t ipadz   = det[3];
  Int_t ipadx   = det[4];

  /*
  // Find out distance d on the plane wrt median phi:
  Float_t d = (ipadx+0.5-kNpadX*0.5)*fgkXPad;

  // The radius r in xy plane:
  //Float_t r = (fgkRmin+fgkRmax)*0.5-0.01+GetHeights(iplate,istrip)+
  //  (ipadz-0.5)*fgkZPad*TMath::Sin(GetAngles(iplate,istrip)/kRaddeg)-0.25; ???
  Float_t r = (fgkRmin+fgkRmax)*0.5-0.01+GetHeights(iplate,istrip)+
    (ipadz-0.5)*fgkZPad*TMath::Sin(GetAngles(iplate,istrip)/kRaddeg);

  // local azimuthal angle in the sector philoc
  Float_t philoc  = TMath::ATan(d/r);
  //if(philoc<0.) philoc = k2PI + philoc;

  // azimuthal angle in the global frame  phi
  Float_t phi   = philoc*kRaddeg+(isector+0.5)*fPhiSec;

  Float_t xCoor = r/TMath::Cos(philoc)*TMath::Cos(phi/kRaddeg);
  */

  // Pad reference frame -> FSTR reference frame
  //  /*
  Float_t posLocal[3] = {0., 0., 0.};
  Float_t step[3] = {-(ipadx+0.5)*fgkXPad, 0., -(ipadz+0.5)*fgkZPad};
  Translation(posLocal,step);

  step[0] = kNpadX*0.5*fgkXPad;
  step[1] = 0.;
  step[2] = kNpadZ*0.5*fgkZPad;
  //  */
  /*
  Float_t posLocal[3] = {(ipadx+0.5)*fgkXPad, 0., (ipadz+0.5)*fgkZPad};
  Float_t step[3]= {kNpadX*0.5*fgkXPad, 0., kNpadZ*0.5*fgkZPad};
  */
  Translation(posLocal,step);

  // FSTR reference frame -> FTOA/B/C = FLTA/B/C reference frame
  Double_t angles[6];
  if      (GetAngles(iplate,istrip) >0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] = GetAngles(iplate,istrip);
    angles[5] = 90.;
  }
  else if (GetAngles(iplate,istrip)==0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.;
    angles[3] = 90.;
    angles[4] =  0;
    angles[5] =  0.;
  }
  else if (GetAngles(iplate,istrip) <0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] =-GetAngles(iplate,istrip);
    angles[5] = 270.;
  }

  InverseRotation(posLocal,angles);

  step[0] = 0.;
  step[1] = -GetHeights(iplate,istrip);
  step[2] =  GetDistances(iplate,istrip);
  Translation(posLocal,step);

  // FTOA = FLTA reference frame -> B071/B074/B075 = BTO1/2/3 reference frame
  angles[0] = 90.;
  angles[1] =  0.;
  angles[2] =  0.;
  angles[3] =  0.;
  angles[4] = 90.;
  angles[5] =270.;

  InverseRotation(posLocal,angles);

  // B071/B074/B075 = BTO1/2/3 reference frame -> ALICE reference frame
  step[0] = 0.;
  step[1] = 0.;
  step[2] = -((fgkRmax+fgkRmin)*0.5);
  Translation(posLocal,step);

  angles[0] = 90.;
  angles[1] = 90.+(isector+0.5)*fPhiSec;
  angles[2] = 0.;
  angles[3] = 0.;
  angles[4] = 90.;
  angles[5] = (isector+0.5)*fPhiSec;

  InverseRotation(posLocal,angles);

  Float_t xCoor = posLocal[0];

  return xCoor;

}
//_____________________________________________________________________________
Float_t AliTOFGeometryV5::GetY(Int_t *det) const
{
  //
  // Returns Y coordinate (cm)
  //

  Int_t isector = det[0];
  Int_t iplate  = det[1];
  Int_t istrip  = det[2];
  Int_t ipadz   = det[3];
  Int_t ipadx   = det[4];

  /*
  // Find out distance d on the plane wrt median phi:
  Float_t d = (ipadx+0.5-kNpadX*0.5)*fgkXPad;

  // The radius r in xy plane:
  //Float_t r = (fgkRmin+fgkRmax)*0.5-0.01+GetHeights(iplate,istrip)+
  //  (ipadz-0.5)*fgkZPad*TMath::Sin(GetAngles(iplate,istrip)/kRaddeg)-0.25; ???
  Float_t r = (fgkRmin+fgkRmax)*0.5-0.01+GetHeights(iplate,istrip)+
    (ipadz-0.5)*fgkZPad*TMath::Sin(GetAngles(iplate,istrip)/kRaddeg);

  // local azimuthal angle in the sector philoc
  Float_t philoc = TMath::ATan(d/r);
  //if(philoc<0.) philoc = k2PI + philoc;

  // azimuthal angle in the global frame  phi
  Float_t phi   = philoc*kRaddeg+(isector+0.5)*fPhiSec;

  Float_t yCoor = r/TMath::Cos(philoc)*TMath::Sin(phi/kRaddeg);
  */

  // Pad reference frame -> FSTR reference frame
  //  /*
  Float_t posLocal[3] = {0., 0., 0.};
  Float_t step[3] = {-(ipadx+0.5)*fgkXPad, 0., -(ipadz+0.5)*fgkZPad};
  Translation(posLocal,step);

  step[0] = kNpadX*0.5*fgkXPad;
  step[1] = 0.;
  step[2] = kNpadZ*0.5*fgkZPad;
  //  */
  /*
  Float_t posLocal[3] = {(ipadx+0.5)*fgkXPad, 0., (ipadz+0.5)*fgkZPad};
  Float_t step[3]= {kNpadX*0.5*fgkXPad, 0., kNpadZ*0.5*fgkZPad};
  */
  Translation(posLocal,step);

  // FSTR reference frame -> FTOA/B/C = FLTA/B/C reference frame

  Double_t angles[6];
  if      (GetAngles(iplate,istrip) >0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] = GetAngles(iplate,istrip);
    angles[5] = 90.;
  }
  else if (GetAngles(iplate,istrip)==0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.;
    angles[3] = 90.;
    angles[4] =  0;
    angles[5] =  0.;
  }
  else if (GetAngles(iplate,istrip) <0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] =-GetAngles(iplate,istrip);
    angles[5] = 270.;
  }

  InverseRotation(posLocal,angles);

  step[0] = 0.;
  step[1] = -GetHeights(iplate,istrip);
  step[2] =  GetDistances(iplate,istrip);
  Translation(posLocal,step);

  // FTOA = FLTA reference frame -> B071/B074/B075 = BTO1/2/3 reference frame
  angles[0] = 90.;
  angles[1] =  0.;
  angles[2] =  0.;
  angles[3] =  0.;
  angles[4] = 90.;
  angles[5] =270.;

  InverseRotation(posLocal,angles);

  // B071/B074/B075 = BTO1/2/3 reference frame -> ALICE reference frame
  step[0] = 0.;
  step[1] = 0.;
  step[2] = -((fgkRmax+fgkRmin)*0.5);
  Translation(posLocal,step);

  angles[0] = 90.;
  angles[1] = 90.+(isector+0.5)*fPhiSec;
  angles[2] = 0.;
  angles[3] = 0.;
  angles[4] = 90.;
  angles[5] = (isector+0.5)*fPhiSec;

  InverseRotation(posLocal,angles);

  Float_t yCoor = posLocal[1];

  return yCoor;

}

//_____________________________________________________________________________
Float_t AliTOFGeometryV5::GetZ(Int_t *det) const
{
  //
  // Returns Z coordinate (cm)
  //

  Int_t isector = det[0];
  Int_t iplate  = det[1];
  Int_t istrip  = det[2];
  Int_t ipadz   = det[3];
  Int_t ipadx   = det[4];

  /*
  Float_t zCoor = GetDistances(iplate,istrip) +
    (0.5-ipadz) * fgkZPad * TMath::Cos(GetAngles(iplate,istrip)*kDegrad);
  */

  // Pad reference frame -> FSTR reference frame
  //  /*
  Float_t posLocal[3] = {0., 0., 0.};
  Float_t step[3] = {-(ipadx+0.5)*fgkXPad, 0., -(ipadz+0.5)*fgkZPad};
  Translation(posLocal,step);

  step[0] = kNpadX*0.5*fgkXPad;
  step[1] = 0.;
  step[2] = kNpadZ*0.5*fgkZPad;
  //  */
  /*
  Float_t posLocal[3] = {(ipadx+0.5)*fgkXPad, 0., (ipadz+0.5)*fgkZPad};
  Float_t step[3]= {kNpadX*0.5*fgkXPad, 0., kNpadZ*0.5*fgkZPad};
  */
  Translation(posLocal,step);

  // FSTR reference frame -> FTOA/B/C = FLTA/B/C reference frame
  Double_t angles[6];
  if      (GetAngles(iplate,istrip) >0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] = GetAngles(iplate,istrip);
    angles[5] = 90.;
  }
  else if (GetAngles(iplate,istrip)==0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.;
    angles[3] = 90.;
    angles[4] =  0;
    angles[5] =  0.;
  }
  else if (GetAngles(iplate,istrip) <0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] =-GetAngles(iplate,istrip);
    angles[5] = 270.;
  }

  InverseRotation(posLocal,angles);

  step[0] = 0.;
  step[1] = -GetHeights(iplate,istrip);
  step[2] =  GetDistances(iplate,istrip);
  Translation(posLocal,step);

  // FTOA = FLTA reference frame -> B071/B074/B075 = BTO1/2/3 reference frame
  angles[0] = 90.;
  angles[1] =  0.;
  angles[2] =  0.;
  angles[3] =  0.;
  angles[4] = 90.;
  angles[5] =270.;

  InverseRotation(posLocal,angles);

  // B071/B074/B075 = BTO1/2/3 reference frame -> ALICE reference frame
  step[0] = 0.;
  step[1] = 0.;
  step[2] = -((fgkRmax+fgkRmin)*0.5);
  Translation(posLocal,step);

  angles[0] = 90.;
  angles[1] = 90.+(isector+0.5)*fPhiSec;
  angles[2] = 0.;
  angles[3] = 0.;
  angles[4] = 90.;
  angles[5] = (isector+0.5)*fPhiSec;

  InverseRotation(posLocal,angles);

  Float_t zCoor = posLocal[2];

  return zCoor;

}

//_____________________________________________________________________________
Int_t AliTOFGeometryV5::GetSector(Float_t *pos) const
{
  //
  // Returns the Sector index 
  //

  //const Float_t khAlWall = 0.1;
  //const Float_t kModuleWallThickness = 0.3;

  Int_t   iSect = -1; 

  Float_t x = pos[0];
  Float_t y = pos[1];
  Float_t z = pos[2];

  Float_t rho = TMath::Sqrt(x*x + y*y);

  //if (!((z>=-fgkMaxhZtof && z<=fgkMaxhZtof) &&
  if (!((z>=-fgkZlenA*0.5 && z<=fgkZlenA*0.5) &&
	(rho>=(fgkRmin) && rho<=(fgkRmax)))) {
    //(rho>=(fgkRmin-0.05)+kModuleWallThickness && rho<=(fgkRmax-0.05)-kModuleWallThickness-khAlWall-kModuleWallThickness))) {
    //AliError("Detector Index could not be determined");
    return iSect;
  }

  Float_t phi = TMath::Pi() + TMath::ATan2(-y,-x);	

  iSect  = (Int_t) (phi*kRaddeg/fPhiSec);
  
  return iSect;

}
//_____________________________________________________________________________

Int_t AliTOFGeometryV5::GetPlate(Float_t *pos) const
{
  //
  // Returns the Plate index 
  //
  const Float_t kInterCentrModBorder1 = 49.5;
  const Float_t kInterCentrModBorder2 = 57.5;
  const Float_t kExterInterModBorder1 = 196.0;
  const Float_t kExterInterModBorder2 = 203.5;

  const Float_t kLengthExInModBorder  = 4.7;
  const Float_t kLengthInCeModBorder  = 7.0;

  //const Float_t khAlWall = 0.1;
  const Float_t kModuleWallThickness = 0.3;
  //const Float_t kHoneycombLayerThickness = 1.5;

  Int_t iPlate=-1;

  Float_t posLocal[3];
  for (Int_t ii=0; ii<3; ii++) posLocal[ii] = pos[ii];

  Int_t isector = GetSector(posLocal);
  if(isector == -1){
    //AliError("Detector Index could not be determined");
    return iPlate;
  }

  // ALICE reference frame -> B071/B074/B075 = BTO1/2/3 reference frame
  Double_t angles[6] = 
    {90., 90.+(isector+0.5)*fPhiSec,
      0., 0.,
     90., (isector+0.5)*fPhiSec
    };
  Rotation(posLocal,angles);

  Float_t step[3] = {0., 0., (fgkRmax+fgkRmin)*0.5};
  Translation(posLocal,step);

  // B071/B074/B075 = BTO1/2/3 reference frame -> FTOA = FLTA reference frame
  angles[0] = 90.;
  angles[1] =  0.;
  angles[2] =  0.;
  angles[3] =  0.;
  angles[4] = 90.;
  angles[5] =270.;

  Rotation(posLocal,angles);

  Float_t yLocal = posLocal[1];
  Float_t zLocal = posLocal[2];

  Float_t deltaRhoLoc  = (fgkRmax-fgkRmin)*0.5 - kModuleWallThickness + yLocal;
  Float_t deltaZetaLoc = TMath::Abs(zLocal);

  Float_t deltaRHOmax = 0.;

  if (TMath::Abs(zLocal)>=kExterInterModBorder1 && TMath::Abs(zLocal)<=kExterInterModBorder2) 
    {
      deltaRhoLoc -= kLengthExInModBorder;
      deltaZetaLoc = kExterInterModBorder2-deltaZetaLoc;
      deltaRHOmax  = (fgkRmax - fgkRmin)*0.5 - kModuleWallThickness - 2.*kLengthExInModBorder; // old 5.35, new 4.8

      if (deltaRhoLoc > deltaZetaLoc*deltaRHOmax/(kInterCentrModBorder2-kInterCentrModBorder1)) {
	if (zLocal<0) iPlate = 0;
	else iPlate = 4;
      }
      else {
	if (zLocal<0) iPlate = 1;
	else iPlate = 3;
      }
    }
  else if (TMath::Abs(zLocal)>=kInterCentrModBorder1 && TMath::Abs(zLocal)<=kInterCentrModBorder2) 
    {
      deltaRhoLoc -= kLengthInCeModBorder;
      deltaZetaLoc = deltaZetaLoc-kInterCentrModBorder1;
      deltaRHOmax = (fgkRmax - fgkRmin)*0.5 - kModuleWallThickness - 2.*kLengthInCeModBorder; // old 0.39, new 0.2

      if (deltaRhoLoc>deltaZetaLoc*deltaRHOmax/(kInterCentrModBorder2-kInterCentrModBorder1)) iPlate = 2;
      else {
	if (zLocal<0) iPlate = 1;
	else iPlate = 3;
      }
    }

  if      (zLocal>-fgkZlenA*0.5/*fgkMaxhZtof*/ && zLocal<-kExterInterModBorder2)       iPlate = 0;
  else if (zLocal>-kExterInterModBorder1       && zLocal<-kInterCentrModBorder2)       iPlate = 1;
  else if (zLocal>-kInterCentrModBorder1       && zLocal< kInterCentrModBorder1)       iPlate = 2;
  else if (zLocal> kInterCentrModBorder2       && zLocal< kExterInterModBorder1)       iPlate = 3;
  else if (zLocal> kExterInterModBorder2       && zLocal< fgkZlenA*0.5/*fgkMaxhZtof*/) iPlate = 4;
  
  return iPlate;

}

//_____________________________________________________________________________
Int_t AliTOFGeometryV5::GetStrip(Float_t *pos) const
{
  //
  // Returns the Strip index 
  //
  const Float_t khhony    = 1.0          ; // heigth of HONY  Layer
  const Float_t khpcby    = 0.08         ; // heigth of PCB   Layer
  const Float_t khrgly    = 0.055        ; // heigth of RED GLASS  Layer
  const Float_t khglfy    = 0.285        ; // heigth of GLASS+FISHLINE  Layer
  const Float_t khcpcby   = 0.16         ; // heigth of PCB  Central Layer
  const Float_t kwcpcbz   = 12.4         ; // z dimension of PCB  Central Layer
  const Float_t khstripy = 2.*khhony+2.*khpcby+4.*khrgly+2.*khglfy+khcpcby;//3.11
  const Float_t kwstripz = kwcpcbz;
  const Float_t klstripx = fgkStripLength;
  
  Int_t iStrip=-1;
   
  Float_t posLocal[3];
  for (Int_t ii=0; ii<3; ii++) posLocal[ii] = pos[ii];
 
  Int_t isector = GetSector(posLocal);
  if(isector == -1){
    //AliError("Detector Index could not be determined");
    return iStrip;}
  Int_t iplate =  GetPlate(posLocal);
  if(iplate == -1){
    //AliError("Detector Index could not be determined");
    return iStrip;} 

  Int_t nstrips=0;
  switch (iplate) {
  case 0:
    nstrips=kNStripC;
    break;
  case 4:
    nstrips=kNStripC;
    break;
  case 1:
    nstrips=kNStripB;
    break;
  case 3:
    nstrips=kNStripB;
    break;
  case 2:
    nstrips=kNStripA;
    break;
  }
  
  // ALICE reference frame -> B071/B074/B075 = BTO1/2/3 reference frame
  Double_t angles[6] = 
    {90., 90.+(isector+0.5)*fPhiSec,
      0., 0.,
     90., (isector+0.5)*fPhiSec
    };
  Rotation(posLocal,angles);

  Float_t step[3] = {0., 0., (fgkRmax+fgkRmin)*0.5};
  Translation(posLocal,step);

  // B071/B074/B075 = BTO1/2/3 reference frame -> FTOA = FLTA reference frame
  angles[0] = 90.;
  angles[1] =  0.;
  angles[2] =  0.;
  angles[3] =  0.;
  angles[4] = 90.;
  angles[5] =270.;

  Rotation(posLocal,angles);

  // FTOA/B/C = FLTA/B/C reference frame -> FSTR reference frame
  Int_t totStrip=0;
  for (Int_t istrip=0; istrip<nstrips; istrip++){

    Float_t posLoc2[3]={posLocal[0],posLocal[1],posLocal[2]};	      

    step[0] = 0.;
    step[1] = GetHeights(iplate,istrip);
    step[2] = -GetDistances(iplate,istrip);
    Translation(posLoc2,step);

    if      (GetAngles(iplate,istrip) >0.) {
      angles[0] = 90.;
      angles[1] =  0.;
      angles[2] = 90.+GetAngles(iplate,istrip);
      angles[3] = 90.;
      angles[4] = GetAngles(iplate,istrip);
      angles[5] = 90.;
    }
    else if (GetAngles(iplate,istrip)==0.) {
      angles[0] = 90.;
      angles[1] =  0.;
      angles[2] = 90.;
      angles[3] = 90.;
      angles[4] =  0;
      angles[5] =  0.;
    }
    else if (GetAngles(iplate,istrip) <0.) {
      angles[0] = 90.;
      angles[1] =  0.;
      angles[2] = 90.+GetAngles(iplate,istrip);
      angles[3] = 90.;
      angles[4] =-GetAngles(iplate,istrip);
      angles[5] = 270.;
    }
    Rotation(posLoc2,angles);

    if ((TMath::Abs(posLoc2[0])<=klstripx*0.5) &&
	(TMath::Abs(posLoc2[1])<=khstripy*0.5) &&
	(TMath::Abs(posLoc2[2])<=kwstripz*0.5)) {
      iStrip = istrip;
      totStrip++;
      for (Int_t jj=0; jj<3; jj++) posLocal[jj]=posLoc2[jj];
      //AliInfo(Form(" posLocal[0] = %f, posLocal[1] = %f, posLocal[2] = %f ", posLocal[0],posLocal[1],posLocal[2]));

      //AliInfo(Form(" GetAngles(%1i,%2i) = %f, pos[0] = %f, pos[1] = %f, pos[2] = %f", iplate, istrip, GetAngles(iplate,istrip), pos[0], pos[1], pos[2]));
      break;
    }

    if (totStrip>1) AliInfo(Form("total strip number found %2i",totStrip));

  }

  return iStrip;
  
}
//_____________________________________________________________________________
Int_t AliTOFGeometryV5::GetPadZ(Float_t *pos) const
{
  //
  // Returns the Pad index along Z 
  //
  //const Float_t klsensmx = kNpadX*fgkXPad;  // length of Sensitive Layer
  //const Float_t khsensmy = 0.05;//0.11;//0.16;// heigth of Sensitive Layer
  //const Float_t kwsensmz = kNpadZ*fgkZPad;  // width of Sensitive Layer

  Int_t iPadZ = -1;

  Float_t posLocal[3];
  for (Int_t ii=0; ii<3; ii++) posLocal[ii] = pos[ii];
 
  Int_t isector = GetSector(posLocal);
  if(isector == -1){
    //AliError("Detector Index could not be determined");
    return iPadZ;}
  Int_t iplate =  GetPlate(posLocal);
  if(iplate == -1){
    //AliError("Detector Index could not be determined");
    return iPadZ;}
  Int_t istrip =  GetStrip(posLocal);
  if(istrip == -1){
    //AliError("Detector Index could not be determined");
    return iPadZ;}

  // ALICE reference frame -> B071/B074/B075 = BTO1/2/3 reference frame
  Double_t angles[6] = 
    {90., 90.+(isector+0.5)*fPhiSec,
      0., 0.,
     90., (isector+0.5)*fPhiSec
    };
  Rotation(posLocal,angles);

  Float_t step[3] = {0., 0., (fgkRmax+fgkRmin)*0.5};
  Translation(posLocal,step);

  // B071/B074/B075 = BTO1/2/3 reference frame -> FTOA = FLTA reference frame
  angles[0] = 90.;
  angles[1] =  0.;
  angles[2] =  0.;
  angles[3] =  0.;
  angles[4] = 90.;
  angles[5] =270.;

  Rotation(posLocal,angles);

  // FTOA/B/C = FLTA/B/C reference frame -> FSTR reference frame
  step[0] = 0.;
  step[1] = GetHeights(iplate,istrip);
  step[2] = -GetDistances(iplate,istrip);
  Translation(posLocal,step);

  if      (GetAngles(iplate,istrip) >0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] = GetAngles(iplate,istrip);
    angles[5] = 90.;
  }
  else if (GetAngles(iplate,istrip)==0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.;
    angles[3] = 90.;
    angles[4] =  0;
    angles[5] =  0.;
  }
  else if (GetAngles(iplate,istrip) <0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] =-GetAngles(iplate,istrip);
    angles[5] = 270.;
  }
  Rotation(posLocal,angles);

  //if (TMath::Abs(posLocal[0])<=klsensmx*0.5 && /*TMath::Abs(posLocal[1])<=khsensmy*0.5+0.005 &&*/ TMath::Abs(posLocal[2])<=kwsensmz*0.5) {
  //if (TMath::Abs(posLocal[1])<=khsensmy*0.5) {

    step[0] =-0.5*kNpadX*fgkXPad;
    step[1] = 0.;
    step[2] =-0.5*kNpadZ*fgkZPad;
    Translation(posLocal,step);

    iPadZ = (Int_t)(posLocal[2]/fgkZPad);
    if (iPadZ==kNpadZ) iPadZ--;
    else if (iPadZ>kNpadZ) iPadZ=-1;

  //}
  // else AliError("Detector Index could not be determined");

  return iPadZ;

}
//_____________________________________________________________________________
Int_t AliTOFGeometryV5::GetPadX(Float_t *pos) const
{
  //
  // Returns the Pad index along X 
  //
  //const Float_t klsensmx = kNpadX*fgkXPad;  // length of Sensitive Layer
  //const Float_t khsensmy = 0.05;//0.11;//0.16;// heigth of Sensitive Layer
  //const Float_t kwsensmz = kNpadZ*fgkZPad;  // width of Sensitive Layer

  Int_t iPadX  = -1;

  Float_t posLocal[3];
  for (Int_t ii=0; ii<3; ii++) posLocal[ii] = pos[ii];
 
  Int_t isector = GetSector(posLocal);
  if(isector == -1){
    //AliError("Detector Index could not be determined");
    return iPadX;}
  Int_t iplate =  GetPlate(posLocal);
  if(iplate == -1){
    //AliError("Detector Index could not be determined");
    return iPadX;} 
  Int_t istrip =  GetStrip(posLocal);
  if(istrip == -1){  
    //AliError("Detector Index could not be determined");
    return iPadX;}

  // ALICE reference frame -> B071/B074/B075 = BTO1/2/3 reference frame
  Double_t angles[6] = 
    {90., 90.+(isector+0.5)*fPhiSec,
      0.,  0.,
     90., (isector+0.5)*fPhiSec
    };
  Rotation(posLocal,angles);

  Float_t step[3] = {0., 0., (fgkRmax+fgkRmin)*0.5};
  Translation(posLocal,step);

  // B071/B074/B075 = BTO1/2/3 reference frame -> FTOA/B/C = FLTA/B/C reference frame
  angles[0] = 90.;
  angles[1] =  0.;
  angles[2] =  0.;
  angles[3] =  0.;
  angles[4] = 90.;
  angles[5] =270.;

  Rotation(posLocal,angles);

  // FTOA/B/C = FLTA/B/C reference frame -> FSTR reference frame
  step[0] = 0.;
  step[1] = GetHeights(iplate,istrip);
  step[2] = -GetDistances(iplate,istrip);
  Translation(posLocal,step);

  if      (GetAngles(iplate,istrip) >0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] = GetAngles(iplate,istrip);
    angles[5] = 90.;
  }
  else if (GetAngles(iplate,istrip)==0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.;
    angles[3] = 90.;
    angles[4] =  0;
    angles[5] =  0.;
  }
  else if (GetAngles(iplate,istrip) <0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] =-GetAngles(iplate,istrip);
    angles[5] = 270.;
  }
  Rotation(posLocal,angles);

  //if (TMath::Abs(posLocal[0])<=klsensmx*0.5 && /*TMath::Abs(posLocal[1])<=khsensmy*0.5+0.005 &&*/ TMath::Abs(posLocal[2])<=kwsensmz*0.5) {
  //if (TMath::Abs(posLocal[1])<=khsensmy*0.5) {

    step[0] =-0.5*kNpadX*fgkXPad;
    step[1] = 0.;
    step[2] =-0.5*kNpadZ*fgkZPad;
    Translation(posLocal,step);

    iPadX = (Int_t)(posLocal[0]/fgkXPad);
    if (iPadX==kNpadX) iPadX--;
    else if (iPadX>kNpadX) iPadX=-1;

  //}
  //else AliError("Detector Index could not be determined");

  return iPadX;

}
//_____________________________________________________________________________

Float_t AliTOFGeometryV5::GetPadDx(Float_t *pos)
{
  //
  // Returns the x coordinate in the Pad reference frame
  //

  Float_t xpad = -2.;

  Float_t posLocal[3];
  for (Int_t ii=0; ii<3; ii++) posLocal[ii] = pos[ii];
 
  Int_t isector = GetSector(posLocal);
  if(isector == -1){
    //AliError("Detector Index could not be determined");
    return xpad;}
  Int_t iplate =  GetPlate(posLocal);
  if(iplate == -1){
    //AliError("Detector Index could not be determined");
    return xpad;} 
  Int_t istrip =  GetStrip(posLocal);
  if(istrip == -1){  
    //AliError("Detector Index could not be determined");
    return xpad;}
  Int_t ipadz =  GetPadZ(posLocal);
  if(ipadz == -1){  
    //AliError("Detector Index could not be determined");
    return xpad;}
  Int_t ipadx =  GetPadX(posLocal);
  if(ipadx == -1){
    //AliError("Detector Index could not be determined");
    return xpad;}

  // ALICE reference frame -> B071/B074/B075 = BTO1/2/3 reference frame
  Double_t angles[6] = 
    {90., 90.+(isector+0.5)*fPhiSec,
      0.,  0.,
     90., (isector+0.5)*fPhiSec
    };
  Rotation(posLocal,angles);

  Float_t step[3] = {0., 0., (fgkRmax+fgkRmin)*0.5};
  Translation(posLocal,step);

  // B071/B074/B075 = BTO1/2/3 reference frame -> FTOA/B/C = FLTA/B/C reference frame
  angles[0] = 90.;
  angles[1] =  0.;
  angles[2] =  0.;
  angles[3] =  0.;
  angles[4] = 90.;
  angles[5] =270.;

  Rotation(posLocal,angles);

  // FTOA/B/C = FLTA/B/C reference frame -> FSTR reference frame
  step[0] = 0.;
  step[1] = GetHeights(iplate,istrip);
  step[2] = -GetDistances(iplate,istrip);
  Translation(posLocal,step);

  if      (GetAngles(iplate,istrip) >0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] = GetAngles(iplate,istrip);
    angles[5] = 90.;
  }
  else if (GetAngles(iplate,istrip)==0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.;
    angles[3] = 90.;
    angles[4] =  0;
    angles[5] =  0.;
  }
  else if (GetAngles(iplate,istrip) <0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] =-GetAngles(iplate,istrip);
    angles[5] = 270.;
  }
  Rotation(posLocal,angles);

  step[0] =-0.5*kNpadX*fgkXPad;
  step[1] = 0.;
  step[2] =-0.5*kNpadZ*fgkZPad;
  Translation(posLocal,step);

  step[0] = (ipadx+0.5)*fgkXPad;
  step[1] = 0.;
  step[2] = (ipadz+0.5)*fgkZPad;
  Translation(posLocal,step);
  
  xpad=posLocal[0];

  return xpad;

}
//_____________________________________________________________________________
Float_t AliTOFGeometryV5::GetPadDy(Float_t *pos)
{
  //
  // Returns the x coordinate in the Pad reference frame
  //

  Float_t ypad = -2.;

  Float_t posLocal[3];
  for (Int_t ii=0; ii<3; ii++) posLocal[ii] = pos[ii];
 
  Int_t isector = GetSector(posLocal);
  if(isector == -1){
    //AliError("Detector Index could not be determined");
    return ypad;}
  Int_t iplate =  GetPlate(posLocal);
  if(iplate == -1){
    //AliError("Detector Index could not be determined");
    return ypad;} 
  Int_t istrip =  GetStrip(posLocal);
  if(istrip == -1){  
    //AliError("Detector Index could not be determined");
    return ypad;}
  Int_t ipadz =  GetPadZ(posLocal);
  if(ipadz == -1){  
    //AliError("Detector Index could not be determined");
    return ypad;}
  Int_t ipadx =  GetPadX(posLocal);
  if(ipadx == -1){
    //AliError("Detector Index could not be determined");
    return ypad;}

  // ALICE reference frame -> B071/B074/B075 = BTO1/2/3 reference frame
  Double_t angles[6] = 
    {90., 90.+(isector+0.5)*fPhiSec,
      0.,  0.,
     90., (isector+0.5)*fPhiSec
    };
  Rotation(posLocal,angles);

  Float_t step[3] = {0., 0., (fgkRmax+fgkRmin)*0.5};
  Translation(posLocal,step);

  // B071/B074/B075 = BTO1/2/3 reference frame -> FTOA/B/C = FLTA/B/C reference frame
  angles[0] = 90.;
  angles[1] =  0.;
  angles[2] =  0.;
  angles[3] =  0.;
  angles[4] = 90.;
  angles[5] =270.;

  Rotation(posLocal,angles);

  // FTOA/B/C = FLTA/B/C reference frame -> FSTR reference frame
  step[0] = 0.;
  step[1] = GetHeights(iplate,istrip);
  step[2] = -GetDistances(iplate,istrip);
  Translation(posLocal,step);

  if      (GetAngles(iplate,istrip) >0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] = GetAngles(iplate,istrip);
    angles[5] = 90.;
  }
  else if (GetAngles(iplate,istrip)==0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.;
    angles[3] = 90.;
    angles[4] =  0;
    angles[5] =  0.;
  }
  else if (GetAngles(iplate,istrip) <0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] =-GetAngles(iplate,istrip);
    angles[5] = 270.;
  }
  Rotation(posLocal,angles);

  step[0] =-0.5*kNpadX*fgkXPad;
  step[1] = 0.;
  step[2] =-0.5*kNpadZ*fgkZPad;
  Translation(posLocal,step);
  
  step[0] = (ipadx+0.5)*fgkXPad;
  step[1] = 0.;
  step[2] = (ipadz+0.5)*fgkZPad;
  Translation(posLocal,step);
  
  ypad=posLocal[1];
  
  return ypad;

}
//_____________________________________________________________________________
Float_t AliTOFGeometryV5::GetPadDz(Float_t *pos)
{
  //
  // Returns the x coordinate in the Pad reference frame
  //

  Float_t zpad = -2.;

  Float_t posLocal[3];
  for (Int_t ii=0; ii<3; ii++) posLocal[ii] = pos[ii];
 
  Int_t isector = GetSector(posLocal);
  if(isector == -1){
    //AliError("Detector Index could not be determined");
    return zpad;}
  Int_t iplate =  GetPlate(posLocal);
  if(iplate == -1){
    //AliError("Detector Index could not be determined");
    return zpad;} 
  Int_t istrip =  GetStrip(posLocal);
  if(istrip == -1){  
    //AliError("Detector Index could not be determined");
    return zpad;}
  Int_t ipadz =  GetPadZ(posLocal);
  if(ipadz == -1){  
    //AliError("Detector Index could not be determined");
    return zpad;}
  Int_t ipadx =  GetPadX(posLocal);
  if(ipadx == -1){
    //AliError("Detector Index could not be determined");
    return zpad;}

  // ALICE reference frame -> B071/B074/B075 = BTO1/2/3 reference frame
  Double_t angles[6] = 
    {90., 90.+(isector+0.5)*fPhiSec,
      0.,  0.,
     90., (isector+0.5)*fPhiSec
    };
  Rotation(posLocal,angles);

  Float_t step[3] = {0., 0., (fgkRmax+fgkRmin)*0.5};
  Translation(posLocal,step);

  // B071/B074/B075 = BTO1/2/3 reference frame -> FTOA/B/C = FLTA/B/C reference frame
  angles[0] = 90.;
  angles[1] =  0.;
  angles[2] =  0.;
  angles[3] =  0.;
  angles[4] = 90.;
  angles[5] =270.;

  Rotation(posLocal,angles);

  // FTOA/B/C = FLTA/B/C reference frame -> FSTR reference frame
  step[0] = 0.;
  step[1] = GetHeights(iplate,istrip);
  step[2] = -GetDistances(iplate,istrip);
  Translation(posLocal,step);

  if      (GetAngles(iplate,istrip) >0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] = GetAngles(iplate,istrip);
    angles[5] = 90.;
  }
  else if (GetAngles(iplate,istrip)==0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.;
    angles[3] = 90.;
    angles[4] =  0;
    angles[5] =  0.;
  }
  else if (GetAngles(iplate,istrip) <0.) {
    angles[0] = 90.;
    angles[1] =  0.;
    angles[2] = 90.+GetAngles(iplate,istrip);
    angles[3] = 90.;
    angles[4] =-GetAngles(iplate,istrip);
    angles[5] = 270.;
  }
  Rotation(posLocal,angles);

  step[0] =-0.5*kNpadX*fgkXPad;
  step[1] = 0.;
  step[2] =-0.5*kNpadZ*fgkZPad;
  Translation(posLocal,step);
  
  step[0] = (ipadx+0.5)*fgkXPad;
  step[1] = 0.;
  step[2] = (ipadz+0.5)*fgkZPad;
  Translation(posLocal,step);

  zpad=posLocal[2];

  return zpad;

}
//_____________________________________________________________________________

void AliTOFGeometryV5::Translation(Float_t *xyz, Float_t translationVector[3]) const
{
  //
  // Return the vector xyz translated by translationVector vector
  //

  Int_t ii=0;

  for (ii=0; ii<3; ii++)
    xyz[ii] -= translationVector[ii];

  return;

}
//_____________________________________________________________________________

void AliTOFGeometryV5::Rotation(Float_t *xyz, Double_t rotationAngles[6]) const
{
  //
  // Return the vector xyz rotated according to the rotationAngles angles
  //

  Int_t ii=0;
  /*
  TRotMatrix *matrix = new TRotMatrix("matrix","matrix", angles[0], angles[1],
				      angles[2], angles[3],
				      angles[4], angles[5]);
  */

  for (ii=0; ii<6; ii++) rotationAngles[ii]*=kDegrad;

  Float_t xyzDummy[3] = {0., 0., 0.};

  for (ii=0; ii<3; ii++) {
    xyzDummy[ii] =
      xyz[0]*TMath::Sin(rotationAngles[2*ii])*TMath::Cos(rotationAngles[2*ii+1]) +
      xyz[1]*TMath::Sin(rotationAngles[2*ii])*TMath::Sin(rotationAngles[2*ii+1]) +
      xyz[2]*TMath::Cos(rotationAngles[2*ii]);
  }

  for (ii=0; ii<3; ii++) xyz[ii]=xyzDummy[ii];

  return;

}
//_____________________________________________________________________________
void AliTOFGeometryV5::InverseRotation(Float_t *xyz, Double_t rotationAngles[6]) const
{
  //
  //
  //

  Int_t ii=0;

  for (ii=0; ii<6; ii++) rotationAngles[ii]*=kDegrad;

  Float_t xyzDummy[3] = {0., 0., 0.};

  xyzDummy[0] =
    xyz[0]*TMath::Sin(rotationAngles[0])*TMath::Cos(rotationAngles[1]) +
    xyz[1]*TMath::Sin(rotationAngles[2])*TMath::Cos(rotationAngles[3]) +
    xyz[2]*TMath::Sin(rotationAngles[4])*TMath::Cos(rotationAngles[5]);
  
  xyzDummy[1] =
    xyz[0]*TMath::Sin(rotationAngles[0])*TMath::Sin(rotationAngles[1]) +
    xyz[1]*TMath::Sin(rotationAngles[2])*TMath::Sin(rotationAngles[3]) +
    xyz[2]*TMath::Sin(rotationAngles[4])*TMath::Sin(rotationAngles[5]);
  
  xyzDummy[2] =
    xyz[0]*TMath::Cos(rotationAngles[0]) +
    xyz[1]*TMath::Cos(rotationAngles[2]) +
    xyz[2]*TMath::Cos(rotationAngles[4]);
  
  for (ii=0; ii<3; ii++) xyz[ii]=xyzDummy[ii];

  return;

}
//_____________________________________________________________________________
void AliTOFGeometryV5::GetVolumePath(Int_t *ind, Char_t *path ) {
  //--------------------------------------------------------------------
  // This function returns the colume path of a given pad 
  //--------------------------------------------------------------------
  Int_t sector = ind[0];
  Char_t  string1[100];
  Char_t  string2[100];
  Char_t  string3[100];
  
  Int_t icopy=-1;
  
  if(sector<3){
    icopy=sector+1;
    sprintf(string1,"/ALIC_1/B077_1/B075_%i/BTO3_1/FTOA_0/FLTA_0",icopy);
  }
  else if(sector<11){
    icopy=sector+3;
    sprintf(string1,"/ALIC_1/B077_1/B071_%i/BTO1_1/FTOA_0/FLTA_0",icopy);
  }
  else if(sector==11 || sector==12){
    icopy=sector-10;
    sprintf(string1,"/ALIC_1/B077_1/B074_%i/BTO2_1/FTOA_0/FLTA_0",icopy);
    if(fHoles)sprintf(string1,"/ALIC_1/B077_1/B074_%i/BTO2_1",icopy);
  }
  else {
    icopy=sector-12;
    sprintf(string1,"/ALIC_1/B077_1/B071_%i/BTO1_1/FTOA_0/FLTA_0",icopy);
  }
  
  Int_t iplate=ind[1];
  Int_t istrip=ind[2];
  if( iplate==0) icopy=istrip; 
  if( iplate==1) icopy=istrip+NStripC(); 
  if( iplate==2) icopy=istrip+NStripC()+NStripB(); 
  if( iplate==3) icopy=istrip+NStripC()+NStripB()+NStripA(); 
  if( iplate==4) icopy=istrip+NStripC()+2*NStripB()+NStripA(); 
  icopy++;
  sprintf(string2,"FSTR_%i",icopy);
  if(fHoles && (sector==11 || sector==12)){
    if(iplate<2)  sprintf(string2,"FTOB_0/FLTB_0/FSTR_%i",icopy);
    if(iplate>2)  sprintf(string2,"FTOC_0/FLTC_0/FSTR_%i",icopy);
  }
 

  Int_t padz = ind[3]+1; 
  Int_t padx = ind[4]+1;
  sprintf(string3,"FPCB_1/FSEN_1/FSEZ_%i/FPAD_%i",padz,padx);
  sprintf(path,"%s/%s/%s",string1,string2,string3); 

}
//_____________________________________________________________________________
void AliTOFGeometryV5::GetPos(Int_t *det, Float_t *pos) 
{
//
// Returns space point coor (x,y,z) (cm)  for Detector 
// Indices  (iSect,iPlate,iStrip,iPadX,iPadZ) 
//
  Char_t path[100];
  GetVolumePath(det,path );
  if (!gGeoManager) {
    printf("ERROR: no TGeo\n");
  }
  gGeoManager->cd(path);
  TGeoHMatrix global;
  global = *gGeoManager->GetCurrentMatrix();
  const Double_t *tr = global.GetTranslation();

  pos[0]=tr[0];  
  pos[1]=tr[1];  
  pos[2]=tr[2];
}
//_____________________________________________________________________________
