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
Revision 1.20.1  2007/05/19 decaro
         Added the following methods:
             GetVolumeIndices(Int_t index, Int_t *det), to get
          the volume indices (sector, plate, strip, padz, padx,
          stored respectively in det[0], det[1], det[2], det[3], det[4])
          from the calibration channel index;
             NStrip(Int_t nPlate), to get the strips number
          per each kind of TOF module.

Revision 1.20  2007/10/08 17:52:55  decaro
hole region in front of PHOS detector: update of sectors' numbers

Revision 1.19  2007/10/04 14:05:09  zampolli
AliTOFGeometryV5 becoming AliTOFGeometry

Revision 1.18  2007/02/19 18:55:26  decaro
Added getter methods for volume path (for Event Display)

Revision 1.17.1  2006/12/15 
         Added method DetToStripRF(...) to get
         a pad corner coordinates in its strip reference frame
         (A.De Caro, M.Di Stefano)
Revision 1.17  2006/08/22 13:30:02  arcelli
removal of effective c++ warnings (C.Zampolli)

Revision 1.16  2006/04/20 22:30:50  hristov
Coding conventions (Annalisa)

Revision 1.15  2006/04/16 22:29:05  hristov
Coding conventions (Annalisa)

Revision 1.14  2006/04/05 08:35:38  hristov
Coding conventions (S.Arcelli, C.Zampolli)

Revision 1.13  2006/03/12 14:37:54  arcelli
 Changes for TOF Reconstruction using TGeo

Revision 1.12  2006/02/28 10:38:00  decaro
AliTOFGeometry::fAngles, AliTOFGeometry::fHeights, AliTOFGeometry::fDistances arrays: dimension definition in the right location

Revision 1.11  2005/12/15 14:17:29  decaro
Correction of some parameter values

Revision 1.10  2005/12/15 08:55:32  decaro
New TOF geometry description (V5) -G. Cara Romeo and A. De Caro

Revision 1.9.1  2005/07/19 A. De Caro
        Created daughter-classes AliTOFGeometryV4 and AliTOFGeometryV5
	=> moved global methods IsInsideThePad, DistanceToPad,
	GetPlate, GetSector, GetStrip, GetPadX, GetPadZ,
	GetX, GetY, GetZ, GetPadDx, GetPadDy and GetPadDz
	in daughter-classes

Revision 1.9  2005/10/20 12:41:35  hristov
Implementation of parallel tracking. It is not the default version, one can use it passing option MI from AliReconstruction to TOF (M.Ivanov)

Revision 1.8  2004/11/29 08:28:01  decaro
Introduction of a new TOF constant (i.e. TDC bin width)

Revision 1.7  2004/11/05 07:20:08  decaro
TOF library splitting and conversion of some printout messages in AliLog schema (T.Kuhr)

Revision 1.6  2004/06/15 15:27:59  decaro
TOF raw data: preliminary implementation and style changes

Revision 1.5  2004/04/20 14:37:22  hristov
Using TMath::Abs instead of fabs, arrays of variable size created/deleted correctly (HP,Sun)

Revision 1.4  2004/04/13 09:42:51  decaro
Track reconstruction code for TOF: updating

Revision 1.3  2003/12/29 18:40:39  hristov
Copy/paste error corrected

Revision 1.2  2003/12/29 17:26:01  hristov
Using enum to initaialize static ints in the header file, the initialization of static floats moved to the implementation file

Revision 1.1  2003/12/29 15:18:03  decaro
TOF geometry updating (addition of AliTOFGeometry)

Revision 0.05  2004/6/11 A.De Caro
        Implement Global method NpadXStrip
        Insert four float constants (originally  in AliTOF class)
Revision 0.04  2004/4/05 S.Arcelli
        Implement Global methods IsInsideThePad 
                                  DistanceToPad 
Revision 0.03  2003/12/14 S.Arcelli
        Set Phi range [-180,180]->[0,360] 
Revision 0.02  2003/12/10 S.Arcelli:
        Implement Global methods GetPos & GetDetID 
Revision 0.01  2003/12/04 S.Arcelli
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TOF Geometry class                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TGeoManager.h"
//#include "TGeoMatrix.h"
#include "TMath.h"

#include "AliConst.h"
#include "AliGeomManager.h"
#include "AliLog.h"

#include "AliTOFGeometry.h"

extern TGeoManager *gGeoManager;

ClassImp(AliTOFGeometry)

const Float_t AliTOFGeometry::fgkZlenA    = 370.6*2.; // length (cm) of the A module
const Float_t AliTOFGeometry::fgkZlenB    = 146.5;    // length (cm) of the B module
const Float_t AliTOFGeometry::fgkZlenC    = 170.45;   // length (cm) of the C module
const Float_t AliTOFGeometry::fgkMaxhZtof = 370.6;    // Max half z-size of TOF (cm)

const Float_t AliTOFGeometry::fgkxTOF     = 372.00;// Inner radius of the TOF for Reconstruction (cm)
const Float_t AliTOFGeometry::fgkRmin     = 371.00;// Inner radius of the TOF (cm)
const Float_t AliTOFGeometry::fgkRmax     = 400.05;// Outer radius of the TOF (cm)

const Float_t AliTOFGeometry::fgkXPad     = 2.5;    // Pad size in the x direction (cm)
const Float_t AliTOFGeometry::fgkZPad     = 3.5;    // Pad size in the z direction (cm)

const Float_t AliTOFGeometry::fgkStripLength = 122.;// Strip Length (rho X phi direction) (cm)

const Float_t AliTOFGeometry::fgkSigmaForTail1= 2.; //Sig1 for simulation of TDC tails 
const Float_t AliTOFGeometry::fgkSigmaForTail2= 0.5;//Sig2 for simulation of TDC tails

const Float_t AliTOFGeometry::fgkPhiSec= 20;//sector Phi width (deg)

const Float_t AliTOFGeometry::fgkTdcBin = 24.4;     // time-of-flight bin width [ps]
const Float_t AliTOFGeometry::fgkToTBin = 48.8;     // time-over-threshold bin width [ps]
const Float_t AliTOFGeometry::fgkBunchCrossingBin = fgkTdcBin * 1024; // bunch-crossing bin width [ps]

const Float_t AliTOFGeometry::fgkSlewTOTMin = 10.; // min TOT for slewing correction [ns]
const Float_t AliTOFGeometry::fgkSlewTOTMax = 16.; // max TOT for slewing correction [ns]

const Float_t AliTOFGeometry::fgkDeadTime = 25E+03;        // Single channel dead time (ps)
const Float_t AliTOFGeometry::fgkMatchingWindow = fgkTdcBin*TMath::Power(2,13); // Matching window  (ps)

const Float_t AliTOFGeometry::fgkAngles[kNPlates][kMaxNstrip] = {
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

/*
const Float_t AliTOFGeometry::fgkHeights[kNPlates][kMaxNstrip] = {
    {-8.2,  -7.5,  -8.2,  -7.7,  -8.1,  -7.6,  -7.7,  -7.7,  -7.7,  -7.7,
     -7.5,  -7.2,  -7.3,  -7.5,  -7.6,  -7.8,  -8.3,  -9.3,  -3.1,   0.0},

    {-7.9,  -8.1,  -8.5,  -9.0, -10.1,  -3.9,  -5.9,  -7.7, -10.1,  -3.6,
     -5.8,  -8.0, -10.4,  -4.4,  -7.2, -10.2,  -4.6,  -7.4, -10.4,   0.0},

    {-2.5, -10.4,  -5.0,  -9.9,  -4.8,  -9.9,  -4.7, -10.2,  -4.7,  -9.9,
     -4.8,  -9.9,  -5.0, -10.4,  -2.5,   0.0,   0.0,   0.0,   0.0,   0.0},

    {-10.4, -7.4,  -4.6, -10.2,  -7.2,  -4.4, -10.4,  -8.0,  -5.8,  -3.6,
     -10.1, -7.7,  -5.9,  -3.9, -10.1,  -9.0,  -8.5,  -8.1,  -7.9,   0.0},

    { -3.1, -9.3,  -8.3,  -7.8,  -7.6,  -7.5,  -7.3,  -7.2,  -7.5,  -7.7,
      -7.7, -7.7,  -7.7,  -7.6,  -8.1,  -7.7,  -8.2,  -7.5,  -8.2,   0.0}
  };
*/
/*
const Float_t AliTOFGeometry::fgkHeights[kNPlates][kMaxNstrip] = {
  {  -8.405, -10.885,  -8.405,  -7.765,  -8.285,  -7.745,  -7.865,  -7.905,  -7.895,  -7.885,
     -7.705,  -7.395,  -7.525,  -7.645, -11.285, -10.355,  -8.365,  -9.385,  -3.255,   0.000 },
  {  -7.905,  -8.235,  -8.605,  -9.045, -10.205,  -3.975,  -5.915,  -7.765, -10.205,  -3.635,
     -5.885,  -8.005, -10.505,  -4.395,  -7.325, -10.235,  -4.655,  -7.495, -10.515,   0.000 },
  {  -2.705, -10.645,  -5.165, -10.095,  -4.995, -10.815,  -4.835, -10.385,  -4.835, -10.815,
     -4.995, -10.095,  -5.165, -10.645,  -2.705,   0.000,   0.000,   0.000,   0.000,   0.000 },
  { -10.515,  -7.495,  -4.655, -10.235,  -7.325,  -4.395, -10.505,  -8.005,  -5.885,  -3.635,
    -10.205,  -7.765,  -5.915,  -3.975, -10.205,  -9.045,  -8.605,  -8.235,  -7.905,   0.000 },
  {  -3.255,  -9.385,  -8.365, -10.355, -11.285,  -7.645,  -7.525,  -7.395,  -7.705,  -7.885,
     -7.895,  -7.905,  -7.865,  -7.745,  -8.285,  -7.765,  -8.405, -10.885,  -8.405,   0.000 }
};
*/


const Float_t AliTOFGeometry::fgkHeights[kNPlates][kMaxNstrip] = {
  { -8.405,  -7.725,  -8.405,  -7.765,  -8.285,  -7.745,  -7.865,  -7.905,  -7.895,  -7.885,
    -7.705,  -7.395,  -7.525,  -7.645,  -7.835,  -7.965,  -8.365,  -9.385,  -3.255,   0.000 },
  { -7.905,  -8.235,  -8.605,  -9.045, -10.205,  -3.975,  -5.915,  -7.765, -10.205,  -3.635,
    -5.885,  -8.005, -10.505,  -4.395,  -7.325, -10.235,  -4.655,  -7.495, -10.515,   0.000 },
  { -2.705, -10.645,  -5.165, -10.095,  -4.995, -10.085,  -4.835, -10.385,  -4.835, -10.085,
    -4.995, -10.095,  -5.165, -10.645,  -2.705,   0.000,   0.000,   0.000,   0.000,   0.000 },
  {-10.515,  -7.495,  -4.655, -10.235,  -7.325,  -4.395, -10.505,  -8.005,  -5.885,  -3.635,
   -10.205,  -7.765,  -5.915,  -3.975, -10.205,  -9.045,  -8.605,  -8.235,  -7.905,   0.000 },
  { -3.255,  -9.385,  -8.365,  -7.965,  -7.835,  -7.645,  -7.525,  -7.395,  -7.705,  -7.885,
    -7.895,  -7.905,  -7.865,  -7.745,  -8.285,  -7.765,  -8.405,  -7.725,  -8.405,   0.000 }
};



const Float_t AliTOFGeometry::fgkDistances[kNPlates][kMaxNstrip] = {
  { 364.14,  354.88,  344.49,  335.31,  325.44,  316.51,  307.11,  297.91,  288.84,  279.89,
    271.20,  262.62,  253.84,  245.20,  236.56,  228.06,  219.46,  210.63,  206.09,    0.00 },
  { 194.57,  186.38,  178.25,  170.13,  161.78,  156.62,  148.10,  139.72,  131.23,  125.87,
    117.61,  109.44,  101.29,   95.46,   87.36,   79.37,   73.17,   65.33,   57.71,    0.00 },
  {  49.28,   41.35,   35.37,   27.91,   21.20,   13.94,    7.06,    0.00,   -7.06,  -13.94,
    -21.20,  -27.91,  -35.37,  -41.35,  -49.28,    0.00,    0.00,    0.00,    0.00,    0.00 },
  { -57.71,  -65.33,  -73.17,  -79.37,  -87.36,  -95.46, -101.29, -109.44, -117.61, -125.87,
   -131.23, -139.72, -148.10, -156.62, -161.78, -170.13, -178.25, -186.38, -194.57,    0.00 },
  {-206.09, -210.63, -219.46, -228.06, -236.56, -245.20, -253.84, -262.62, -271.20, -279.89,
   -288.84, -297.91, -307.11, -316.51, -325.44, -335.31, -344.49, -354.88, -364.14,    0.00 }
};

/*
const Float_t AliTOFGeometry::fgkDistances[kNPlates][kMaxNstrip] = {
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
*/
//_____________________________________________________________________________
AliTOFGeometry::AliTOFGeometry():
  fHoles(1)
{
  //
  // AliTOFGeometry default constructor
  //

}

//_____________________________________________________________________________
AliTOFGeometry::~AliTOFGeometry()
{
  //
  // AliTOFGeometry destructor
  //
}
//_____________________________________________________________________________
void AliTOFGeometry::ImportGeometry(){
  TGeoManager::Import("geometry.root");
}
//_____________________________________________________________________________
void AliTOFGeometry::GetPosPar(Int_t *det, Float_t *pos) const
{
//
// Returns space point coor (x,y,z) (cm)  for Detector 
// Indices  (iSect,iPlate,iStrip,iPadX,iPadZ) 
//

  pos[0]=GetX(det);  
  pos[1]=GetY(det);  
  pos[2]=GetZ(det);
  
}
//_____________________________________________________________________________
void AliTOFGeometry::GetDetID( Float_t *pos, Int_t *det) const
{
 //
 // Returns Detector Indices (iSect,iPlate,iStrip,iPadX,iPadZ) 
 // space point coor (x,y,z) (cm)  


  det[0]=GetSector(pos);  
  det[1]=GetPlate(pos);  
  det[2]=GetStrip(pos);
  det[3]=GetPadZ(pos);
  det[4]=GetPadX(pos);
  
}
//_____________________________________________________________________________

void AliTOFGeometry::DetToStripRF(Int_t nPadX, Int_t nPadZ, Float_t &x,  Float_t &z) const
{
  //
  // Returns the local coordinates (x, z) in strip reference frame
  // for the bottom corner of the pad number (nPadX, nPadZ)
  //
  /*
  const Float_t xCenterStrip = kNpadX * fgkXPad / 2.;
  const Float_t zCenterStrip = kNpadZ * fgkZPad / 2.;

  const Float_t xCenterPad = nPadX*fgkXPad + fgkXPad / 2.;
  const Float_t zCenterPad = nPadZ*fgkZPad + fgkZPad / 2.;

  x = xCenterPad - xCenterStrip;
  z = zCenterPad - zCenterStrip;
  */


  x = (nPadX - kNpadX*0.5) * fgkXPad;
  z = (nPadZ - kNpadZ*0.5) * fgkZPad;


}
//_____________________________________________________________________________
Float_t AliTOFGeometry::DistanceToPadPar(Int_t *det, const Float_t * const pos, Float_t *dist3d) const
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
  Float_t angle   = phi*kRaddeg-( Int_t (kRaddeg*phi/fgkPhiSec) + 0.5)*fgkPhiSec;
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
  Float_t padAngle = padPhi*kRaddeg-( Int_t (padPhi*kRaddeg/fgkPhiSec)+ 0.5) * fgkPhiSec;
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
Bool_t AliTOFGeometry::IsInsideThePadPar(Int_t *det, const Float_t * const pos) const
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

  const Float_t kPadDepth = 0.5;//0.05;//0.11;//0.16;//          // heigth of Sensitive Layer

  //Transform pos into Sector Frame

  Float_t x = pos[0];
  Float_t y = pos[1];
  Float_t z = pos[2];

  Float_t radius = TMath::Sqrt(x*x+y*y);
  Float_t phi = TMath::Pi()+TMath::ATan2(-y,-x);	

  //  Get the local angle in the sector philoc
  Float_t angle = phi*kRaddeg-( Int_t (kRaddeg*phi/fgkPhiSec) + 0.5) *fgkPhiSec;
  Float_t xs = radius*TMath::Cos(angle/kRaddeg);
  Float_t ys = radius*TMath::Sin(angle/kRaddeg);
  Float_t zs = z;

  // Do the same for the selected pad

  Float_t g[3];
  GetPosPar(det,g);

  Float_t padRadius = TMath::Sqrt(g[0]*g[0]+g[1]*g[1]);
  Float_t padPhi = TMath::Pi()+TMath::ATan2(-g[1],-g[0]);	

  //  Get the local angle in the sector philoc
  Float_t padAngle = padPhi*kRaddeg-( Int_t (padPhi*kRaddeg/fgkPhiSec)+ 0.5) * fgkPhiSec; 
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

  if(TMath::Abs(xr)<=kPadDepth*0.5 && TMath::Abs(yr)<= (fgkXPad*0.5) && TMath::Abs(zr)<= (fgkZPad*0.5))
    isInside=true;
  return isInside;

}
//_____________________________________________________________________________
Bool_t AliTOFGeometry::IsInsideThePad(TGeoHMatrix mat, const Float_t * const pos, Float_t *dist3d) const
{
  //
  // Returns true if space point with coor pos (x,y,z) [cm] falls inside
  // pad identified by the matrix mat. In case dist3d!=0, dist3d vector
  // has been filled with the 3D distance between the impact point on
  // the pad and the pad centre (in the reference frame of the TOF pad
  // identified by the matrix mat).
  //

  const Float_t kPadDepth = 0.5;      // heigth of Sensitive Layer

  Double_t posg[3];
  posg[0] = pos[0];
  posg[1] = pos[1];
  posg[2] = pos[2];

  // from ALICE global reference system
  // towards TOF pad reference system
  Double_t posl[3] = {0., 0., 0.};
  mat.MasterToLocal(posg,posl);

  Float_t xr = posl[0];
  Float_t yr = posl[1];
  Float_t zr = posl[2];

  Bool_t isInside = false;
  if (TMath::Abs(yr)<= kPadDepth*0.5 &&
      TMath::Abs(xr)<= fgkXPad*0.5 &&
      TMath::Abs(zr)<= fgkZPad*0.5)
    isInside = true;

  if (dist3d) {
    //Double_t padl[3] = {0., 0., 0.};
    dist3d[0] = posl[0]/* - padl[0]*/;
    dist3d[1] = posl[1]/* - padl[1]*/;
    dist3d[2] = posl[2]/* - padl[2]*/;

    /*
    Double_t padg[3] = {0., 0., 0.};
    // from TOF pad local reference system
    // towards ALICE global reference system
    TGeoHMatrix inverse = mat.Inverse();
    inverse.MasterToLocal(padl,padg);

    // returns the 3d distance
    // between the impact point on the pad
    // and the pad centre (in the ALICE global reference frame)
    dist3d[0] = posg[0] - padg[0];
    dist3d[1] = posg[1] - padg[1];
    dist3d[2] = posg[2] - padg[2];
    */
  }
 
  return isInside;

}
//_____________________________________________________________________________
void AliTOFGeometry::GetVolumePath(const Int_t * const ind, Char_t *path ) {
  //--------------------------------------------------------------------
  // This function returns the colume path of a given pad 
  //--------------------------------------------------------------------
  Int_t sector = ind[0];
  Char_t  string1[100];
  Char_t  string2[100];
  Char_t  string3[100];
  
  Int_t icopy=-1;
  icopy=sector;
 
  sprintf(string1,"/ALIC_1/B077_1/BSEGMO%i_1/BTOF%i_1",icopy,icopy);
  
  Int_t iplate=ind[1];
  Int_t istrip=ind[2];
  if( iplate==0) icopy=istrip; 
  if( iplate==1) icopy=istrip+NStripC(); 
  if( iplate==2) icopy=istrip+NStripC()+NStripB(); 
  if( iplate==3) icopy=istrip+NStripC()+NStripB()+NStripA(); 
  if( iplate==4) icopy=istrip+NStripC()+2*NStripB()+NStripA(); 
  icopy++;
  sprintf(string2,"FTOA_0/FLTA_0/FSTR_%i",icopy);
  if(fHoles && (sector==13 || sector==14 || sector==15)){
    if(iplate<2)  sprintf(string2,"FTOB_0/FLTB_0/FSTR_%i",icopy);
    if(iplate>2)  sprintf(string2,"FTOC_0/FLTC_0/FSTR_%i",icopy);
  }
 
  Int_t padz = ind[3]+1; 
  Int_t padx = ind[4]+1;
  sprintf(string3,"FPCB_1/FSEN_1/FSEZ_%i/FPAD_%i",padz,padx);
  sprintf(path,"%s/%s/%s",string1,string2,string3); 

}
//_____________________________________________________________________________
void AliTOFGeometry::GetVolumePath(Int_t sector, Char_t *path ){
  //--------------------------------------------------------------------
  // This function returns the colume path of a given sector 
  //--------------------------------------------------------------------

  Char_t string[100];

  Int_t icopy = sector;

  sprintf(string,"/ALIC_1/B077_1/BSEGMO%i_1/BTOF%i_1",icopy,icopy);
  sprintf(path,"%s",string);

}
//_____________________________________________________________________________
void AliTOFGeometry::GetVolumePath(Int_t sector, Int_t plate, Int_t strip, Char_t *path ) {
  //--------------------------------------------------------------------
  // This function returns the colume path of a given strip 
  //--------------------------------------------------------------------

  Char_t string1[100];
  Char_t string2[100];
  Char_t string3[100];
  
  Int_t icopy = sector;

  sprintf(string1,"/ALIC_1/B077_1/BSEGMO%i_1/BTOF%i_1",icopy,icopy);
  
  if(plate==0) icopy=strip; 
  if(plate==1) icopy=strip+NStripC(); 
  if(plate==2) icopy=strip+NStripC()+NStripB(); 
  if(plate==3) icopy=strip+NStripC()+NStripB()+NStripA(); 
  if(plate==4) icopy=strip+NStripC()+2*NStripB()+NStripA(); 
  icopy++;
  sprintf(string2,"FTOA_0/FLTA_0/FSTR_%i",icopy);
  if(fHoles && (sector==13 || sector==14 || sector==15)){
    if(plate<2)  sprintf(string2,"FTOB_0/FLTB_0/FSTR_%i",icopy);
    if(plate>2)  sprintf(string2,"FTOC_0/FLTC_0/FSTR_%i",icopy);
  }

  sprintf(string3,"FPCB_1/FSEN_1");
  sprintf(path,"%s/%s/%s",string1,string2,string3); 

}
//_____________________________________________________________________________
void AliTOFGeometry::GetPos(Int_t *det, Float_t *pos) 
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
Int_t AliTOFGeometry::GetPlate(const Float_t * const pos) const
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
    {90., 90.+(isector+0.5)*fgkPhiSec,
      0., 0.,
     90., (isector+0.5)*fgkPhiSec
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

  if      (zLocal>-fgkZlenA*0.5          && zLocal<-kExterInterModBorder2) iPlate = 0;
  else if (zLocal>-kExterInterModBorder1 && zLocal<-kInterCentrModBorder2) iPlate = 1;
  else if (zLocal>-kInterCentrModBorder1 && zLocal< kInterCentrModBorder1) iPlate = 2;
  else if (zLocal> kInterCentrModBorder2 && zLocal< kExterInterModBorder1) iPlate = 3;
  else if (zLocal> kExterInterModBorder2 && zLocal< fgkZlenA*0.5)          iPlate = 4;
  
  return iPlate;

}

//_____________________________________________________________________________
Int_t AliTOFGeometry::GetSector(const Float_t * const pos) const
{
  //
  // Returns the Sector index 
  //

  Int_t   iSect = -1; 

  Float_t x = pos[0];
  Float_t y = pos[1];
  Float_t z = pos[2];

  Float_t rho = TMath::Sqrt(x*x + y*y);

  if (!((z>=-fgkZlenA*0.5 && z<=fgkZlenA*0.5) &&
	(rho>=(fgkRmin) && rho<=(fgkRmax)))) {
    //AliError("Detector Index could not be determined");
    return iSect;
  }

  Float_t phi = TMath::Pi() + TMath::ATan2(-y,-x);	

  iSect  = (Int_t) (phi*kRaddeg/fgkPhiSec);
  
  return iSect;

}
//_____________________________________________________________________________
Int_t AliTOFGeometry::GetStrip(const Float_t * const pos) const
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
  AliDebug(1,Form("  posLocal[0] = %f, posLocal[1] = %f, posLocal[2] = %f ",
		  posLocal[0],posLocal[1],posLocal[2]));

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
  case 4:
    nstrips=kNStripC;
    break;
  case 1:
  case 3:
    nstrips=kNStripB;
    break;
  case 2:
    nstrips=kNStripA;
    break;
  }
  
  // ALICE reference frame -> B071/B074/B075 = BTO1/2/3 reference frame
  Double_t angles[6] = 
    {90., 90.+(isector+0.5)*fgkPhiSec,
      0., 0.,
     90., (isector+0.5)*fgkPhiSec
    };
  Rotation(posLocal,angles);
  AliDebug(1,Form("  posLocal[0] = %f, posLocal[1] = %f, posLocal[2] = %f ",
		  posLocal[0],posLocal[1],posLocal[2]));

  Float_t step[3] = {0., 0., (fgkRmax+fgkRmin)*0.5};
  Translation(posLocal,step);
  AliDebug(1,Form("  posLocal[0] = %f, posLocal[1] = %f, posLocal[2] = %f ",
		  posLocal[0],posLocal[1],posLocal[2]));

  // B071/B074/B075 = BTO1/2/3 reference frame -> FTOA = FLTA reference frame
  angles[0] = 90.;
  angles[1] =  0.;
  angles[2] =  0.;
  angles[3] =  0.;
  angles[4] = 90.;
  angles[5] =270.;

  Rotation(posLocal,angles);
  AliDebug(1,Form("  posLocal[0] = %f, posLocal[1] = %f, posLocal[2] = %f ",
		  posLocal[0],posLocal[1],posLocal[2]));

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
    AliDebug(1,Form(" strip %2d:  posLoc2[0] = %f, posLoc2[1] = %f, posLoc2[2] = %f ",
		    istrip, posLoc2[0],posLoc2[1],posLoc2[2]));

    if ((TMath::Abs(posLoc2[0])<=klstripx*0.5) &&
	(TMath::Abs(posLoc2[1])<=khstripy*0.5) &&
	(TMath::Abs(posLoc2[2])<=kwstripz*0.5)) {
      iStrip = istrip;
      totStrip++;
      for (Int_t jj=0; jj<3; jj++) posLocal[jj]=posLoc2[jj];
      AliDebug(2,Form(" posLocal[0] = %f, posLocal[1] = %f, posLocal[2] = %f ",
		      posLocal[0],posLocal[1],posLocal[2]));

      AliDebug(2,Form(" GetAngles(%1i,%2i) = %f, pos[0] = %f, pos[1] = %f, pos[2] = %f",
		      iplate, istrip, GetAngles(iplate,istrip), pos[0], pos[1], pos[2]));
      break;
    }

    if (totStrip>1) AliInfo(Form("total strip number found %2i",totStrip));

  }

  return iStrip;
  
}
//_____________________________________________________________________________
Int_t AliTOFGeometry::GetPadZ(const Float_t * const pos) const
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
    {90., 90.+(isector+0.5)*fgkPhiSec,
      0., 0.,
     90., (isector+0.5)*fgkPhiSec
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

  step[0] =-0.5*kNpadX*fgkXPad;
  step[1] = 0.;
  step[2] =-0.5*kNpadZ*fgkZPad;
  Translation(posLocal,step);

  iPadZ = (Int_t)(posLocal[2]/fgkZPad);
  if (iPadZ==kNpadZ) iPadZ--;
  else if (iPadZ>kNpadZ) iPadZ=-1;

  return iPadZ;

}
//_____________________________________________________________________________
Int_t AliTOFGeometry::GetPadX(const Float_t * const pos) const
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
    {90., 90.+(isector+0.5)*fgkPhiSec,
      0.,  0.,
     90., (isector+0.5)*fgkPhiSec
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

  iPadX = (Int_t)(posLocal[0]/fgkXPad);
  if (iPadX==kNpadX) iPadX--;
  else if (iPadX>kNpadX) iPadX=-1;

  return iPadX;

}
//_____________________________________________________________________________
Float_t AliTOFGeometry::GetX(const Int_t * const det) const
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
  Float_t phi   = philoc*kRaddeg+(isector+0.5)*fgkPhiSec;

  Float_t xCoor = r/TMath::Cos(philoc)*TMath::Cos(phi/kRaddeg);
  */

  // Pad reference frame -> FSTR reference frame
  Float_t posLocal[3] = {0., 0., 0.};
  Float_t step[3] = {-(ipadx+0.5)*fgkXPad, 0., -(ipadz+0.5)*fgkZPad};
  Translation(posLocal,step);

  step[0] = kNpadX*0.5*fgkXPad;
  step[1] = 0.;
  step[2] = kNpadZ*0.5*fgkZPad;
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
  angles[1] = 90.+(isector+0.5)*fgkPhiSec;
  angles[2] = 0.;
  angles[3] = 0.;
  angles[4] = 90.;
  angles[5] = (isector+0.5)*fgkPhiSec;

  InverseRotation(posLocal,angles);

  Float_t xCoor = posLocal[0];

  return xCoor;

}
//_____________________________________________________________________________
Float_t AliTOFGeometry::GetY(const Int_t * const det) const
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
  Float_t phi   = philoc*kRaddeg+(isector+0.5)*fgkPhiSec;

  Float_t yCoor = r/TMath::Cos(philoc)*TMath::Sin(phi/kRaddeg);
  */

  // Pad reference frame -> FSTR reference frame
  Float_t posLocal[3] = {0., 0., 0.};
  Float_t step[3] = {-(ipadx+0.5)*fgkXPad, 0., -(ipadz+0.5)*fgkZPad};
  Translation(posLocal,step);

  step[0] = kNpadX*0.5*fgkXPad;
  step[1] = 0.;
  step[2] = kNpadZ*0.5*fgkZPad;
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
  angles[1] = 90.+(isector+0.5)*fgkPhiSec;
  angles[2] = 0.;
  angles[3] = 0.;
  angles[4] = 90.;
  angles[5] = (isector+0.5)*fgkPhiSec;

  InverseRotation(posLocal,angles);

  Float_t yCoor = posLocal[1];

  return yCoor;

}

//_____________________________________________________________________________
Float_t AliTOFGeometry::GetZ(const Int_t * const det) const
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
  Float_t posLocal[3] = {0., 0., 0.};
  Float_t step[3] = {-(ipadx+0.5)*fgkXPad, 0., -(ipadz+0.5)*fgkZPad};
  Translation(posLocal,step);

  step[0] = kNpadX*0.5*fgkXPad;
  step[1] = 0.;
  step[2] = kNpadZ*0.5*fgkZPad;
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
  angles[1] = 90.+(isector+0.5)*fgkPhiSec;
  angles[2] = 0.;
  angles[3] = 0.;
  angles[4] = 90.;
  angles[5] = (isector+0.5)*fgkPhiSec;

  InverseRotation(posLocal,angles);

  Float_t zCoor = posLocal[2];

  return zCoor;

}
//_____________________________________________________________________________

void AliTOFGeometry::DetToSectorRF(Int_t vol[5], Double_t **coord)
{
  //
  // Returns the local coordinates (x, y, z) in sector reference frame
  // for the 4 corners of each sector pad (vol[1], vol[2], vol[3], vol[4])
  //

  if (!gGeoManager) printf("ERROR: no TGeo\n");

  // ALICE -> TOF Sector
  Char_t path1[100]="";
  GetVolumePath(vol[0],path1);
  gGeoManager->cd(path1);
  TGeoHMatrix aliceToSector;
  aliceToSector = *gGeoManager->GetCurrentMatrix();

  // TOF Sector -> ALICE
  //TGeoHMatrix sectorToALICE = aliceToSector.Inverse();

  // ALICE -> TOF Pad
  Char_t path2[100]="";
  GetVolumePath(vol,path2);
  gGeoManager->cd(path2);
  TGeoHMatrix aliceToPad;
  aliceToPad = *gGeoManager->GetCurrentMatrix();

  // TOF Pad -> ALICE
  TGeoHMatrix padToALICE = aliceToPad.Inverse();

  // TOF Pad -> TOF Sector
  TGeoHMatrix padToSector = padToALICE*aliceToSector;

  // TOF Sector -> TOF Pad
  //TGeoHMatrix sectorToPad = sectorToALICE*aliceToPad;

  // coordinates of the pad bottom corner
  Double_t **cornerPad = new Double_t*[4];
  for (Int_t ii=0; ii<4; ii++) cornerPad[ii] = new Double_t[3];

  cornerPad[0][0] = -fgkXPad/2.;
  cornerPad[0][1] =  0.;
  cornerPad[0][2] = -fgkZPad/2.;

  cornerPad[1][0] =  fgkXPad/2.;
  cornerPad[1][1] =  0.;
  cornerPad[1][2] = -fgkZPad/2.;

  cornerPad[2][0] =  fgkXPad/2.;
  cornerPad[2][1] =  0.;
  cornerPad[2][2] =  fgkZPad/2.;

  cornerPad[3][0] = -fgkXPad/2.;
  cornerPad[3][1] =  0.;
  cornerPad[3][2] =  fgkZPad/2.;

  for(Int_t aa=0; aa<4; aa++) for(Int_t bb=0; bb<3; bb++) coord[aa][bb]=0.;

  for (Int_t jj=0; jj<4; jj++) padToSector.MasterToLocal(&cornerPad[jj][0], &coord[jj][0]);

  delete cornerPad;

  //sectorToPad.LocalToMaster(cornerPad, coord);

}
//_____________________________________________________________________________
Float_t AliTOFGeometry::GetPadDx(const Float_t * const pos)
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
    {90., 90.+(isector+0.5)*fgkPhiSec,
      0.,  0.,
     90., (isector+0.5)*fgkPhiSec
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
Float_t AliTOFGeometry::GetPadDy(const Float_t * const pos)
{
  //
  // Returns the y coordinate in the Pad reference frame
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
    {90., 90.+(isector+0.5)*fgkPhiSec,
      0.,  0.,
     90., (isector+0.5)*fgkPhiSec
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
Float_t AliTOFGeometry::GetPadDz(const Float_t * const pos)
{
  //
  // Returns the z coordinate in the Pad reference frame
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
    {90., 90.+(isector+0.5)*fgkPhiSec,
      0.,  0.,
     90., (isector+0.5)*fgkPhiSec
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

void AliTOFGeometry::Translation(Float_t *xyz, Float_t translationVector[3]) const
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

void AliTOFGeometry::Rotation(Float_t *xyz, Double_t rotationAngles[6]) const
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
void AliTOFGeometry::InverseRotation(Float_t *xyz, Double_t rotationAngles[6]) const
{
  //
  // Rotates the vector xyz acordint to the rotationAngles
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

Int_t AliTOFGeometry::GetIndex(const Int_t * const detId)
{
  //Retrieve calibration channel index 
  Int_t isector = detId[0];
  if (isector >= kNSectors){
    printf("Wrong sector number in TOF (%d) !\n",isector);
    return -1;
  }
  Int_t iplate = detId[1];
  if (iplate >= kNPlates){
    printf("Wrong plate number in TOF (%d) !\n",iplate);
    return -1;
  }
  Int_t istrip = detId[2];
  Int_t stripOffset = GetStripNumberPerSM(iplate,istrip);
  if (stripOffset==-1) {
    printf("Wrong strip number per SM in TOF (%d) !\n",stripOffset);
    return -1;
  }

  Int_t ipadz = detId[3];
  Int_t ipadx = detId[4];

  Int_t idet = ((2*(kNStripC+kNStripB)+kNStripA)*kNpadZ*kNpadX)*isector +
               (stripOffset*kNpadZ*kNpadX)+
	       (kNpadX)*ipadz+
	        ipadx;
  return idet;
}
//_____________________________________________________________________________

void AliTOFGeometry::GetVolumeIndices(Int_t index, Int_t *detId)
{
  //
  // Retrieve volume indices from the calibration channel index 
  //

  detId[0] = index/NpadXStrip()/NStripXSector();

  Int_t dummyStripPerModule = 
    ( index - ( NStripXSector()*NpadXStrip()*detId[0]) ) / NpadXStrip();
  if (dummyStripPerModule<kNStripC) {
    detId[1] = 0;
    detId[2] = dummyStripPerModule;
  }
  else if (dummyStripPerModule>=kNStripC && dummyStripPerModule<kNStripC+kNStripB) {
    detId[1] = 1;
    detId[2] = dummyStripPerModule-kNStripC;
  }
  else if (dummyStripPerModule>=kNStripC+kNStripB && dummyStripPerModule<kNStripC+kNStripB+kNStripA) {
    detId[1] = 2;
    detId[2] = dummyStripPerModule-kNStripC-kNStripB;
  }
  else if (dummyStripPerModule>=kNStripC+kNStripB+kNStripA && dummyStripPerModule<kNStripC+kNStripB+kNStripA+kNStripB) {
    detId[1] = 3;
    detId[2] = dummyStripPerModule-kNStripC-kNStripB-kNStripA;
  }
  else if (dummyStripPerModule>=kNStripC+kNStripB+kNStripA+kNStripB && dummyStripPerModule<NStripXSector()) {
    detId[1] = 4;
    detId[2] = dummyStripPerModule-kNStripC-kNStripB-kNStripA-kNStripB;
  }

  Int_t padPerStrip = ( index - ( NStripXSector()*NpadXStrip()*detId[0]) ) - dummyStripPerModule*NpadXStrip();

  detId[3] = padPerStrip / kNpadX; // padZ
  detId[4] = padPerStrip - detId[3]*kNpadX; // padX

}
//_____________________________________________________________________________

Int_t AliTOFGeometry::NStrip(Int_t nPlate)
{
  //
  // Returns the strips number for the plate number 'nPlate'
  //

  Int_t nStrips = kNStripC;

  switch(nPlate) {
  case 2:
    nStrips = kNStripA;
    break;
  case 1:
  case 3:
    nStrips = kNStripB;
    break;
  case 0:
  case 4:
  default:
    nStrips = kNStripC;
    break;
  }

  return nStrips;

}
//-------------------------------------------------------------------------

UShort_t AliTOFGeometry::GetAliSensVolIndex(Int_t isector, Int_t iplate, Int_t istrip) const
{
  //
  // Get the index of the TOF alignable volume in the AliGeomManager order.
  //

  Int_t index = GetStripNumber(isector, iplate, istrip);

  UShort_t volIndex = AliGeomManager::LayerToVolUID(AliGeomManager::kTOF,index);

  return volIndex;

}
//-------------------------------------------------------------------------

Int_t AliTOFGeometry::GetStripNumber(Int_t isector, Int_t iplate, Int_t istrip)
{
  //
  // Get the serial number of the TOF strip number istrip [0,14/18],
  //   in the module number iplate [0,4],
  //   in the TOF SM number isector [0,17].
  // This number will range in [0,1637].
  //

  Bool_t check = (isector >= kNSectors);

  if (check)
    printf("E-AliTOFGeometry::GetStripNumber: Wrong sector number in TOF (%d)!\n",isector);

  Int_t index = -1;
  Int_t stripInSM = GetStripNumberPerSM(iplate, istrip);
  if (!check && stripInSM!=-1)
    index = (2*(kNStripC+kNStripB)+kNStripA)*isector + stripInSM;

  return index;

}
//-------------------------------------------------------------------------

void AliTOFGeometry::GetStripAndModule(Int_t iStripPerSM, Int_t &iplate, Int_t &istrip)
{
  //
  // Convert the serial number of the TOF strip number iStripPerSM [0,90]
  // in module number iplate [0,4] and strip number istrip [0,14/18].
  //

  if (iStripPerSM<0 || iStripPerSM>=kNStripC+kNStripB+kNStripA+kNStripB+kNStripC) {
    iplate = -1;
    istrip = -1;
  }
  else if (iStripPerSM<kNStripC) {
    iplate = 0;
    istrip = iStripPerSM;
  }
  else if (iStripPerSM>=kNStripC && iStripPerSM<kNStripC+kNStripB) {
    iplate = 1;
    istrip = iStripPerSM-kNStripC;
  }
  else if (iStripPerSM>=kNStripC+kNStripB && iStripPerSM<kNStripC+kNStripB+kNStripA) {
    iplate = 2;
    istrip = iStripPerSM-kNStripC-kNStripB;
  }
  else if (iStripPerSM>=kNStripC+kNStripB+kNStripA && iStripPerSM<kNStripC+kNStripB+kNStripA+kNStripB) {
    iplate = 3;
    istrip = iStripPerSM-kNStripC-kNStripB-kNStripA;
  }
  else if (iStripPerSM>=kNStripC+kNStripB+kNStripA+kNStripB && iStripPerSM<kNStripC+kNStripB+kNStripA+kNStripB+kNStripC) {
    iplate = 4;
    istrip = iStripPerSM-kNStripC-kNStripB-kNStripA-kNStripB;
  }


}
//-------------------------------------------------------------------------

Int_t AliTOFGeometry::GetStripNumberPerSM(Int_t iplate, Int_t istrip)
{
  //
  // Get the serial number of the TOF strip number istrip [0,14/18],
  //   in the module number iplate [0,4].
  // This number will range in [0,90].
  //

  Int_t index = -1;

  Bool_t check = (
		  (iplate<0 || iplate>=kNPlates)
		  ||
		  (
		   (iplate==2 && (istrip<0 || istrip>=kNStripA))
		   ||
		   (iplate!=2 && (istrip<0 || istrip>=kNStripC))
		   )
		  );

  if (iplate<0 || iplate>=kNPlates)
    printf("E-AliTOFGeometry::GetStripNumberPerSM: Wrong plate number in TOF (%1d)!\n",iplate);

  if (
      (iplate==2 && (istrip<0 || istrip>=kNStripA))
      ||
      (iplate!=2 && (istrip<0 || istrip>=kNStripC))
      )
    printf("E-AliTOFGeometry::GetStripNumberPerSM: Wrong strip number in TOF "
	   "(strip=%2d in the plate=%1d)!\n",istrip,iplate);

  Int_t stripOffset = 0;
  switch (iplate) {
  case 0:
    stripOffset = 0;
    break;
  case 1:
    stripOffset = kNStripC;
    break;
  case 2:
    stripOffset = kNStripC+kNStripB;
    break;
  case 3:
    stripOffset = kNStripC+kNStripB+kNStripA;
    break;
  case 4:
    stripOffset = kNStripC+kNStripB+kNStripA+kNStripB;
    break;
  };

  if (!check) index = stripOffset + istrip;

  return index;

}
//-------------------------------------------------------------------------

void AliTOFGeometry::PadRF2TrackingRF(Float_t *ctrackPos, Float_t *differenceT)
{
  //
  // To convert the 3D distance ctrackPos, referred to the ALICE RF,
  // into the 3D distance differenceT, referred to the tracking RF
  // in case ctrakPos belongs to a TOF sensitive volume.
  //

  for (Int_t ii=0; ii<3; ii++) differenceT[ii] = 999.;

  AliDebug(1,Form(" track position in ALICE global Ref. frame -> %f, %f, %f",
		  ctrackPos[0],ctrackPos[1],ctrackPos[2]));

  Int_t detId[5] = {-1,-1,-1,-1,-1};

  detId[0] = GetSector(ctrackPos);
  if (detId[0]==-1) {
    AliWarning(Form("This point does not belong to any TOF sector"));
    return;
  }

  detId[1] = GetPlate(ctrackPos);
  if (detId[1]==-1) {
    AliWarning(Form("This point does not belong to any TOF module"));
    return;
  }

  detId[2] = GetStrip(ctrackPos);
  if (detId[2]==-1) {
    AliWarning(Form("This point does not belong to any TOF strip"));
    return;
  }

  detId[3] = GetPadZ(ctrackPos);
  if (detId[3]==-1) {
    AliWarning(Form("This point does not belong to any TOF pad-row"));
    return;
  }

  detId[4] = GetPadX(ctrackPos);
  if (detId[4]==-1) {
    AliWarning(Form("This point does not belong to any TOF pad"));
    return;
  }


  UShort_t alignableStripIndex =
    GetAliSensVolIndex(detId[0],detId[1],detId[2]);
  AliDebug(1,Form(" sector = %2d, plate = %1d, strip = %2d (padZ = %1d, padX = %2d) "
		  "---> stripIndex = %4d",
		  detId[0], detId[1], detId[2], detId[3], detId[4], alignableStripIndex));

  // pad centre coordinates in the strip ref. frame
  Double_t padCentreL[3] = {(detId[4]-AliTOFGeometry::NpadX()/2)*AliTOFGeometry::XPad()
			    +AliTOFGeometry::XPad()/2.,
			    0.,
			    (detId[3]-AliTOFGeometry::NpadZ()/2)*AliTOFGeometry::XPad()
			    +AliTOFGeometry::XPad()/2.};
  // pad centre coordinates in the strip tracking frame
  Double_t padCentreT[3] = {0., 0., 0.};
  TGeoHMatrix l2t = *AliGeomManager::GetTracking2LocalMatrix(alignableStripIndex);
  l2t.MasterToLocal(padCentreL,padCentreT);


  Char_t path[100];
  // pad centre coordinates in its ref. frame
  Double_t padCentreL2[3] = {0., 0., 0.};
  // pad centre coordinates in the ALICE global ref. frame
  Double_t padCentreG[3] = {0., 0., 0.};
  GetVolumePath(detId,path);
  gGeoManager->cd(path);
  TGeoHMatrix g2l = *gGeoManager->GetCurrentMatrix();
  TGeoHMatrix l2g = g2l.Inverse();
  l2g.MasterToLocal(padCentreL2,padCentreG);


  Char_t path2[100];
  // strip centre coordinates in its ref. frame
  Double_t stripCentreL[3] = {0., 0., 0.};
  // strip centre coordinates in the ALICE global ref. frame
  Double_t stripCentreG[3] = {0., 0., 0.};
  GetVolumePath(detId[0],detId[1],detId[2],path2);
  gGeoManager->cd(path2);
  TGeoHMatrix g2lb = *gGeoManager->GetCurrentMatrix();
  TGeoHMatrix l2gb = g2lb.Inverse();
  l2gb.MasterToLocal(stripCentreL,stripCentreG);

  TGeoHMatrix g2t = 0;
  AliGeomManager::GetTrackingMatrix(alignableStripIndex, g2t);

  // track position in the ALICE global ref. frame
  Double_t posG[3];
  for (Int_t ii=0; ii<3; ii++) posG[ii] = (Double_t)ctrackPos[ii];

  // strip centre coordinates in the tracking ref. frame
  Double_t stripCentreT[3] = {0., 0., 0.};
  // track position in the tracking ref. frame
  Double_t posT[3] = {0., 0., 0.};
  g2t.MasterToLocal(posG,posT);
  g2t.MasterToLocal(stripCentreG,stripCentreT);

  for (Int_t ii=0; ii<3; ii++)
    AliDebug(1,Form(" track position in ALICE global and tracking RFs -> posG[%d] = %f --- posT[%d] = %f",
		    ii, posG[ii], ii, posT[ii]));
  for (Int_t ii=0; ii<3; ii++)
    AliDebug(1,Form(" pad centre coordinates in its, the ALICE global and tracking RFs -> "
		    "padCentreL[%d] = %f --- padCentreG[%d] = %f --- padCentreT[%d] = %f",
		    ii, padCentreL[ii],
		    ii, padCentreG[ii],
		    ii, padCentreT[ii]));
  for (Int_t ii=0; ii<3; ii++)
    AliDebug(1,Form(" strip centre coordinates in its, the ALICE global and tracking RFs -> "
		    "stripCentreL[%d] = %f --- stripCentreG[%d] = %f --- stripCentreT[%d] = %f",
		    ii, stripCentreL[ii],
		    ii, stripCentreG[ii],
		    ii, stripCentreT[ii]));
  for (Int_t ii=0; ii<3; ii++)
    AliDebug(1,Form(" difference between the track position and the pad centre in the tracking RF "
		    "-> posT[%d]-padCentreT[%d] = %f",
		    ii,ii,
		    posT[ii]-padCentreT[ii]));

  for (Int_t ii=0; ii<3; ii++) differenceT[ii] = (Float_t)(posT[ii]-padCentreT[ii]);

}
//-------------------------------------------------------------------------

Int_t AliTOFGeometry::GetTOFsupermodule(const Int_t index)
{
  // Return the TOF supermodule where TOF channel index is located

  if (index<0 || index>=NPadXSector()*NSectors()) return -1;
  else return index/NpadXStrip()/NStripXSector();

}
