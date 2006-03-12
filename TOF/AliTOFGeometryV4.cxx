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
Revision 1.2  2006/02/28 10:38:00  decaro
AliTOFGeometry::fAngles, AliTOFGeometry::fHeights, AliTOFGeometry::fDistances arrays: dimension definition in the right location

Revision 1.1  2005/12/15 08:55:33  decaro
New TOF geometry description (V5) -G. Cara Romeo and A. De Caro

Revision 0.1  2005/07/19 A. De Caro
        Modify Global methods IsInsideThePad & DistanceToPad
               according to the PPR TOF geometry
        Implement Global  methods GetPadDx & GetPadDy & GetPadDz
        Modify Global methods GetDetID & GetPlate & GetSector &
                              GetStrip & GetPadX & GetPadZ
               according to the PPR TOF geometry
        Modify Global methods GetPos & GetX & GetY & GetZ
               according to the PPR TOF geometry
*/

#include <stdlib.h>
#include <Riostream.h>
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TOF Geometry class (PPR version)                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliConst.h"

#include "AliTOFGeometry.h"
#include "AliTOFGeometryV4.h"

ClassImp(AliTOFGeometryV4)

const Int_t AliTOFGeometryV4::kNStripC      = 20;       // number of strips in C type module

const Float_t AliTOFGeometryV4::fgkZlenA    = 106.0;    // length (cm) of the A module
const Float_t AliTOFGeometryV4::fgkZlenB    = 141.0;    // length (cm) of the B module
const Float_t AliTOFGeometryV4::fgkZlenC    = 177.5;    // length (cm) of the C module
const Float_t AliTOFGeometryV4::fgkMaxhZtof = 371.5;    // Max half z-size of TOF (cm)

const Float_t AliTOFGeometryV4::fgkDeadBndX = 1.0;      // Dead Boundaries of a Strip along X direction (length) (cm)
const Float_t AliTOFGeometryV4::fgkDeadBndZ = 1.5;      // Dead Boundaries of a Strip along Z direction (width) (cm)
const Float_t AliTOFGeometryV4::fgkOverSpc = 15.3;      // Space available for sensitive layers in radial direction (cm)

const Float_t AliTOFGeometryV4::fgkDprecMin = 0.0000075;//num.prec.tolerance on Thmin 
const Float_t AliTOFGeometryV4::fgkDprecMax = 0.0000100;//num.prec.tolerance on Thma 
const Float_t AliTOFGeometryV4::fgkDprecCen = 0.0000005;//num.prec.tolerance on <Theta> 

const Float_t AliTOFGeometryV4::fgkxTOF     = 371.;     // Inner radius of the TOF for Reconstruction (cm)
const Float_t AliTOFGeometryV4::fgkRmin     = 370.;     // Inner radius of the TOF (cm)
const Float_t AliTOFGeometryV4::fgkRmax     = 399.;     // Outer radius of the TOF (cm)

//_____________________________________________________________________________
AliTOFGeometryV4::AliTOFGeometryV4()
  :AliTOFGeometry()
{
  //
  // AliTOFGeometryV4 default constructor
  //

  AliTOFGeometry::kNStripC   = kNStripC;         // number of strips in C type module

  AliTOFGeometry::kZlenA    = fgkZlenA;          // length (cm) of the A module
  AliTOFGeometry::kZlenB    = fgkZlenB;          // length (cm) of the B module
  AliTOFGeometry::kZlenC    = fgkZlenC;          // length (cm) of the C module
  AliTOFGeometry::kMaxhZtof = fgkMaxhZtof;       // Max half z-size of TOF (cm)

  AliTOFGeometry::fgkxTOF   = fgkxTOF;           // Inner radius of the TOF for Reconstruction (cm)
  AliTOFGeometry::fgkRmin   = fgkRmin;           // Inner radius of the TOF (cm)
  AliTOFGeometry::fgkRmax   = fgkRmax;           // Outer radius of the TOF (cm)

  Init();

}

//_____________________________________________________________________________
AliTOFGeometryV4::~AliTOFGeometryV4()
{
  //
  // AliTOFGeometryV4 destructor
  //

}
//_____________________________________________________________________________
void AliTOFGeometryV4::ImportGeometry(){
  TGeoManager::Import("geometry.root");
}
//_____________________________________________________________________________
void AliTOFGeometryV4::Init()
{
  //
  // Initialize strip Tilt Angles and Heights
  //
  // Strips Tilt Angles
 
  Float_t const kangles[kNPlates][kMaxNstrip] ={

 {44.494, 43.725, 42.946, 42.156, 41.357, 40.548, 39.729, 38.899, 
  38.060, 37.211, 36.353, 35.484, 34.606, 33.719, 32.822, 31.916, 
  31.001, 30.077, 29.144, 28.202 },

 {26.884, 25.922, 24.952, 23.975, 22.989, 22.320, 21.016, 20.309,
  19.015, 18.270, 16.989, 16.205, 14.941, 14.117, 12.871, 12.008,
  10.784, 9.8807, 8.681, 0.0 },

 { 7.5835, 6.4124, 5.4058, 4.2809, 3.2448,  2.1424, 1.078, -0., -1.078, 
  -2.1424, -3.2448, -4.2809, -5.4058, -6.4124, -7.5835, 0.0, 0.0, 0.0,
  0.0, 0.0 },
  
 {-8.681, -9.8807, -10.784, -12.008, -12.871, -14.117, -14.941, -16.205,
  -16.989, -18.27, -19.015, -20.309, -21.016, -22.32, -22.989,
   -23.975, -24.952, -25.922, -26.884, 0. },
  
 {-28.202, -29.144, -30.077, -31.001, -31.916, -32.822, -33.719, -34.606,
  -35.484, -36.353, -37.211, -38.06, -38.899, -39.729, -40.548,
   -41.357, -42.156, -42.946, -43.725, -44.494 }};


  //Strips Heights

   Float_t const kheights[kNPlates][kMaxNstrip]= {

  {-5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5,
   -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5 },
  
  {-6.3, -7.1, -7.9, -8.7, -9.5, -3, -9.5,   -3, -9.5,   -3, 
   -9.5, -3.0, -9.5, -3.0, -9.5, -3, -9.5,   -3,   -9 , 0.},
  
  {  -3,   -9, -4.5,   -9, -4.5,     -9, -4.5,   -9, -4.5,   -9, 
     -4.5,   -9, -4.5,   -9,   -3,   0.0, 0.0, 0.0, 0.0, 0.0 },
  
  {  -9,   -3, -9.5,   -3, -9.5, -3, -9.5,   -3, -9.5,   -3, -9.5,
     -3, -9.5,   -3, -9.5,  -8.7, -7.9, -7.1, -6.3, 0. },
  
  {-5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5,
   -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5, -5.5 }};

   // Deposit in fAngles, fHeights

   for (Int_t iplate = 0; iplate < kNPlates; iplate++) {
     for (Int_t istrip = 0; istrip < kMaxNstrip; istrip++) {
       AliTOFGeometry::fAngles[iplate][istrip]   = kangles[iplate][istrip];
       AliTOFGeometry::fHeights[iplate][istrip]  = kheights[iplate][istrip];
     }
   }

}

//_____________________________________________________________________________
Float_t AliTOFGeometryV4::DistanceToPadPar(Int_t *det, Float_t *pos, Float_t *dist3d) 
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
  Float_t phi=TMath::ATan2(y,x);	
  if(phi<0) phi=2.*TMath::Pi()+phi;
  //  Get the local angle in the sector philoc
  Float_t angle   = phi*kRaddeg-( Int_t (kRaddeg*phi/fPhiSec) + 0.5)*fPhiSec;
  Float_t xs = radius*TMath::Cos(angle/kRaddeg);
  Float_t ys = radius*TMath::Sin(angle/kRaddeg);
  Float_t zs = z;

  // Do the same for the selected pad

  Float_t g[3];
  GetPosPar(det,g);

  Float_t padRadius = TMath::Sqrt(g[0]*g[0]+g[1]*g[1]);
  Float_t padPhi=TMath::ATan2(g[1],g[0]);	
  if(padPhi<0) padPhi=2.*TMath::Pi()+padPhi;
  //  Get the local angle in the sector philoc
  Float_t padAngle   = padPhi*kRaddeg-( Int_t (padPhi*kRaddeg/fPhiSec)+ 0.5) * fPhiSec; 
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
Bool_t AliTOFGeometryV4::IsInsideThePadPar(Int_t *det, Float_t *pos) 
{
//
// Returns true if space point with coor pos (x,y,z) (cm) falls 
// inside pad with Detector Indices idet (iSect,iPlate,iStrip,iPadX,iPadZ) 
//

  Bool_t isInside=false; 

  //Transform pos into Sector Frame

  Float_t x = pos[0];
  Float_t y = pos[1];
  Float_t z = pos[2];

  Float_t radius = TMath::Sqrt(x*x+y*y);
  Float_t phi=TMath::ATan2(y,x);	
  if(phi<0) phi=2.*TMath::Pi()+phi;
  //  Get the local angle in the sector philoc
  Float_t angle   = phi*kRaddeg-( Int_t (kRaddeg*phi/fPhiSec) + 0.5) *fPhiSec;
  Float_t xs = radius*TMath::Cos(angle/kRaddeg);
  Float_t ys = radius*TMath::Sin(angle/kRaddeg);
  Float_t zs = z;

  // Do the same for the selected pad

  Float_t g[3];
  GetPosPar(det,g);

  Float_t padRadius = TMath::Sqrt(g[0]*g[0]+g[1]*g[1]);
  Float_t padPhi=TMath::ATan2(g[1],g[0]);	
  if(padPhi<0) padPhi=2.*TMath::Pi()+padPhi;
  //  Get the local angle in the sector philoc
  Float_t padAngle   = padPhi*kRaddeg-( Int_t (padPhi*kRaddeg/fPhiSec)+ 0.5) * fPhiSec; 
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

  if(TMath::Abs(xr)<=0.75 && TMath::Abs(yr)<= (fgkXPad*0.5) && TMath::Abs(zr)<= (fgkZPad*0.5))
    isInside=true; 
  return isInside;

}


//_____________________________________________________________________________
Float_t AliTOFGeometryV4::DistanceToPad(Int_t *det, TGeoHMatrix mat, Float_t *pos, Float_t *dist3d) 
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
  vecl[2]=-veclr[2];
  //Take into account reflections
  if(det[1]>2){
    vecl[1]=-veclr[0];
    vecl[2]= veclr[2];
  }	

  Float_t dist = TMath::Sqrt(vecl[0]*vecl[0]+vecl[1]*vecl[1]+vecl[2]*vecl[2]);


  if (dist3d){
    dist3d[0] = vecl[0];
    dist3d[1] = vecl[1];
    dist3d[2] = vecl[2];
  }

  return dist;

}


//_____________________________________________________________________________
Bool_t AliTOFGeometryV4::IsInsideThePad( Int_t *det, TGeoHMatrix mat, Float_t *pos) 
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
  mat.MasterToLocal(vecg,veclr);  
  vecl[0]=veclr[1];
  vecl[1]=veclr[0];
  vecl[2]=-veclr[2];
  //Take into account reflections
  if(det[1]>2){
    vecl[1]=-veclr[0];
    vecl[2]= veclr[2];
  }	

  Float_t xr = vecl[0];
  Float_t yr = vecl[1];
  Float_t zr = vecl[2];

  Bool_t isInside=false; 
  if(TMath::Abs(xr)<= khsensmy*0.5 && TMath::Abs(yr)<= (fgkXPad*0.5) && TMath::Abs(zr)<= (fgkZPad*0.5))
    isInside=true; 
  return isInside;

}
//_____________________________________________________________________________
Float_t AliTOFGeometryV4::GetX(Int_t *det) 
{
  //
  // Returns X coordinate (cm)
  //

  Int_t isector = det[0];
  Int_t iplate  = det[1];
  Int_t istrip  = det[2];
  Int_t ipadz   = det[3];
  Int_t ipadx   = det[4];

  // Find out distance d on the plane wrt median phi:
  Float_t d = (ipadx+0.5)*fgkXPad-(kNpadX*fgkXPad)*0.5;

  // The radius r in xy plane:
  Float_t r = (fgkRmin+fgkRmax)/2.+fHeights[iplate][istrip]+
    (ipadz-0.5)*fgkZPad*TMath::Sin(fAngles[iplate][istrip]/kRaddeg)-0.25;

  // local azimuthal angle in the sector philoc
  Float_t philoc   = TMath:: ATan(d/r);

  // azimuthal angle in the global frame  phi
  Float_t phi      = philoc*kRaddeg+(isector+0.5 )*fPhiSec;                    

  Float_t xCoor    = r/TMath::Cos(philoc)*TMath::Cos(phi/kRaddeg);

  return xCoor;

}
//_____________________________________________________________________________
Float_t AliTOFGeometryV4::GetY(Int_t *det) 
{
  //
  // Returns Y coordinate (cm)
  //

  Int_t isector = det[0];
  Int_t iplate  = det[1];
  Int_t istrip  = det[2];
  Int_t ipadz   = det[3];
  Int_t ipadx   = det[4];

  // Find out distance d on the plane wrt median phi:
  Float_t d = (ipadx+0.5)*fgkXPad-(kNpadX*fgkXPad)*0.5;

  // The radius r in xy plane:
  Float_t r = (fgkRmin+fgkRmax)/2.+fHeights[iplate][istrip]+
    (ipadz-0.5)*fgkZPad*TMath::Sin(fAngles[iplate][istrip]/kRaddeg)-0.25;

  // local azimuthal angle in the sector philoc
  Float_t philoc   = TMath:: ATan(d/r);

  // azimuthal angle in the global frame  phi
  Float_t phi      = philoc*kRaddeg+(isector+0.5 )*fPhiSec;                    

  Float_t yCoor    = r/TMath::Cos(philoc)*TMath::Sin(phi/kRaddeg);

  return yCoor;

}

//_____________________________________________________________________________
Float_t AliTOFGeometryV4::GetZ(Int_t *det) 
{
  //
  // Returns Z coordinate (cm)
  //

  Int_t iplate  = det[1];
  Int_t istrip  = det[2];
  Int_t ipadz   = det[3];

  // The radius r in xy plane:
  Float_t r = (fgkRmin+fgkRmax)/2.+fHeights[iplate][istrip];

  Float_t zCoor = r*TMath::Tan(0.5*TMath::Pi()-GetStripTheta(iplate,istrip))-
         (ipadz-0.5)*fgkZPad*TMath::Cos(fAngles[iplate][istrip]/kRaddeg);
  return zCoor;

}

//_____________________________________________________________________________
Int_t AliTOFGeometryV4::GetSector(Float_t *pos) 
{
  //
  // Returns the Sector index 
  //

  Int_t   iSect = -1; 

  Float_t x = pos[0];
  Float_t y = pos[1];

  Float_t phi     =  TMath::ATan2(y,x);	
  if(phi<0.) phi=2.*TMath::Pi()+phi;
  iSect  = (Int_t) (phi*kRaddeg/fPhiSec);

  return iSect;

}

//_____________________________________________________________________________
Int_t AliTOFGeometryV4::GetPadX(Float_t *pos) 
{
  //
  // Returns the Pad index along X 
  //

  Int_t iPadX  = -1;

  Float_t x = pos[0];
  Float_t y = pos[1];
  Float_t z = pos[2];

  Int_t isector = GetSector(pos);
  if(isector == -1){  
    AliError("Detector Index could not be determined");
    return iPadX;}
  Int_t iplate =  GetPlate(pos);
  if(iplate == -1){  
    AliError("Detector Index could not be determined");
    return iPadX;} 
  Int_t istrip =  GetStrip(pos);
  if(istrip == -1){  
    AliError("Detector Index could not be determined");
    return iPadX;}


  Float_t rho=TMath::Sqrt(x*x+y*y);
  Float_t phi =  TMath::ATan2(y,x);	
  if(phi<0.) phi=2.*TMath::Pi()+phi;
 
  // Get the local angle in the sector philoc
  Float_t philoc   = phi*kRaddeg-(isector+0.5)*fPhiSec;
  philoc*=TMath::Pi()/180.;
  // theta projected on the median of the sector
  Float_t theta = TMath::ATan2(rho*TMath::Cos(philoc),z);
  // The radius r in xy plane:
  Float_t r   = (fgkRmin+fgkRmax)/2.+fHeights[iplate][istrip]+
               (theta-GetStripTheta(iplate, istrip))/
    (GetMaxStripTheta(iplate, istrip)-GetMinStripTheta(iplate, istrip))
   * 2.*fgkZPad*TMath::Sin(fAngles[iplate][istrip]/kRaddeg)-0.25;

  // Find out distance projected onto the strip plane 
  Float_t d = (r*TMath::Tan(philoc)+(kNpadX*fgkXPad)*0.5);

  iPadX  =  (Int_t) ( d/fgkXPad);  
  return iPadX;

}
//_____________________________________________________________________________
Int_t AliTOFGeometryV4::GetPlate(Float_t *pos) 
{
  //
  // Returns the Plate index 
  //
  Int_t iPlate=-1;

  Int_t isector = GetSector(pos);
  if(isector == -1){  
    AliError("Detector Index could not be determined");
    return iPlate;}
 
  Float_t x = pos[0];
  Float_t y = pos[1];
  Float_t z = pos[2];

  Float_t rho=TMath::Sqrt(x*x+y*y);
  Float_t phi=TMath::ATan2(y,x);	
  if(phi<0) phi=2.*TMath::Pi()+phi;
  // Get the local angle in the sector philoc
  Float_t philoc   = phi*kRaddeg-(isector+0.5)*fPhiSec;
  philoc*=TMath::Pi()/180.;
  // theta projected on the median of the sector
  Float_t theta=TMath::ATan2(rho*TMath::Cos(philoc),z);

  for (Int_t i=0; i<kNPlates; i++){
    if ( GetMaxPlateTheta(i) >= theta && 
         GetMinPlateTheta(i) <= theta)iPlate=i;
  }

  return iPlate;

}

//_____________________________________________________________________________
Int_t AliTOFGeometryV4::GetStrip(Float_t *pos) 
{
  //
  // Returns the Strip index 
  //

  Int_t iStrip=-1;


  Int_t isector = GetSector(pos);
  if(isector == -1){  
    AliError("Detector Index could not be determined");
    return iStrip;}
  Int_t iplate =  GetPlate(pos);
  if(iplate == -1){  
    AliError("Detector Index could not be determined");
    return iStrip;} 


  Float_t x = pos[0];
  Float_t y = pos[1];
  Float_t z = pos[2];

  Int_t nstrips=0;
  if(iplate==0 || iplate == 4)nstrips=kNStripC;
  if(iplate==1 || iplate == 3)nstrips=kNStripB;
  if(iplate==2)               nstrips=kNStripA;

  Float_t rho=TMath::Sqrt(x*x+y*y);
  Float_t phi=TMath::ATan2(y,x);	
  if(phi<0) phi=2.*TMath::Pi()+phi;
  // Get the local angle in the sector philoc
  Float_t philoc   = phi*kRaddeg-(isector+0.5)*fPhiSec;
  philoc*=TMath::Pi()/180.;
  // theta projected on the median of the sector
  Float_t theta=TMath::ATan2(rho*TMath::Cos(philoc),z);

  for (Int_t istrip=0; istrip<nstrips; istrip++){

    if( 
       GetMaxStripTheta(iplate,istrip) >= theta 
       &&  
       GetMinStripTheta(iplate,istrip) <= theta ) iStrip = istrip;
   
  }

  return iStrip;

}
//_____________________________________________________________________________
Int_t AliTOFGeometryV4::GetPadZ(Float_t *pos) 
{
  //
  // Returns the Pad index along Z 
  //
  Int_t iPadZ = -1;

  Int_t isector = GetSector(pos);
  if(isector == -1){  
    AliError("Detector Index could not be determined");
    return iPadZ;}
  Int_t iplate =  GetPlate(pos);
  if(iplate == -1){  
    AliError("Detector Index could not be determined");
    return iPadZ;} 
  Int_t istrip =  GetStrip(pos);
  if(istrip == -1){  
    AliError("Detector Index could not be determined");
    return iPadZ;}


  Float_t x = pos[0];
  Float_t y = pos[1];
  Float_t z = pos[2];

  Float_t rho=TMath::Sqrt(x*x+y*y);
  Float_t phi=TMath::ATan2(y,x);	
  if(phi<0) phi=2.*TMath::Pi()+phi;
  Float_t philoc   = phi*kRaddeg-(isector+0.5)*fPhiSec;
  philoc*=TMath::Pi()/180.;
  Float_t theta=TMath::ATan2(rho*TMath::Cos(philoc),z);

  if (theta >= GetStripTheta(iplate, istrip))iPadZ=1;
  else iPadZ=0;

  return iPadZ;

}
//_____________________________________________________________________________
Float_t AliTOFGeometryV4::GetMinPlateTheta(Int_t iPlate) 
{
  //
  // Returns the minimum theta angle of a given plate iPlate (rad)
  //
  

  Int_t index=0;

  Float_t delta =0.;
  if(iPlate==0)delta = -1. ;
  if(iPlate==1)delta = -0.5;
  if(iPlate==3)delta = +0.5;
  if(iPlate==4)delta = +1. ;

  Float_t z=(fgkRmin+2.)*TMath::Tan(fAngles[iPlate][index]/kRaddeg)+delta;
  Float_t r=(fgkRmin+fgkRmax)/2.+fHeights[iPlate][index];
  z =z+fgkZPad*TMath::Cos(fAngles[iPlate][index]/kRaddeg);
  r =r-fgkZPad*TMath::Sin(fAngles[iPlate][index]/kRaddeg);

  Float_t thmin = 0.5*TMath::Pi()-TMath::ATan(z/r)-fgkDprecMin;
  return thmin;

}
//_____________________________________________________________________________
Float_t AliTOFGeometryV4::GetMaxPlateTheta(Int_t iPlate) 
{
  //
  // Returns the maximum theta angle of a given plate iPlate (rad)
  
  Int_t index=0;
  if(iPlate==0 ||iPlate == 4)index=kNStripC-1;
  if(iPlate==1 ||iPlate == 3)index=kNStripB-1;
  if(iPlate==2)              index=kNStripA-1;

  Float_t delta =0.;
  if(iPlate==0)delta = -1. ;
  if(iPlate==1)delta = -0.5;
  if(iPlate==3)delta = +0.5;
  if(iPlate==4)delta = +1. ;

  Float_t z=(fgkRmin+2.)*TMath::Tan(fAngles[iPlate][index]/kRaddeg)+delta;
  Float_t r=(fgkRmin+fgkRmax)/2.+fHeights[iPlate][index];
  z =z-fgkZPad*TMath::Cos(fAngles[iPlate][index]/kRaddeg);
  r= r+fgkZPad*TMath::Sin(fAngles[iPlate][index]/kRaddeg);

  Float_t thmax    = 0.5*TMath::Pi()-TMath::ATan(z/r)+fgkDprecMax;

  return thmax;

}
//_____________________________________________________________________________
Float_t  AliTOFGeometryV4::GetMaxStripTheta(Int_t iPlate, Int_t iStrip) 
{
  //
  // Returns the maximum theta angle of a given strip iStrip (rad)
  //


  Float_t delta =0.;
  if(iPlate==0)delta = -1. ;
  if(iPlate==1)delta = -0.5;
  if(iPlate==3)delta = +0.5;
  if(iPlate==4)delta = +1. ;

  Float_t r =(fgkRmin+fgkRmax)/2.+fHeights[iPlate][iStrip];
  Float_t z =(fgkRmin+2.)*TMath::Tan(fAngles[iPlate][iStrip]/kRaddeg)+delta;
  z = z-fgkZPad*TMath::Cos(fAngles[iPlate][iStrip]/kRaddeg);
  r = r+fgkZPad*TMath::Sin(fAngles[iPlate][iStrip]/kRaddeg);
  Float_t thmax =0.5*TMath::Pi()-TMath::ATan(z/r)+fgkDprecMax;
  return thmax;

}
//_____________________________________________________________________________
Float_t  AliTOFGeometryV4::GetMinStripTheta(Int_t iPlate, Int_t iStrip) 
{
  //
  // Returns the minimum theta angle of a given Strip iStrip (rad)
  //
  

  Float_t delta =0.;
  if(iPlate==0)delta = -1. ;
  if(iPlate==1)delta = -0.5;
  if(iPlate==3)delta = +0.5;
  if(iPlate==4)delta = +1. ;


  Float_t r =(fgkRmin+fgkRmax)/2.+fHeights[iPlate][iStrip];
  Float_t z =(fgkRmin+2.)*TMath::Tan(fAngles[iPlate][iStrip]/kRaddeg)+delta;
  z =z+fgkZPad*TMath::Cos(fAngles[iPlate][iStrip]/kRaddeg);
  r =r-fgkZPad*TMath::Sin(fAngles[iPlate][iStrip]/kRaddeg);
  Float_t thmin =0.5*TMath::Pi()-TMath::ATan(z/r)-fgkDprecMin;

  return thmin;

}
//_____________________________________________________________________________
Float_t  AliTOFGeometryV4::GetStripTheta(Int_t iPlate, Int_t iStrip) 
{
  //
  // returns the median theta angle of a given strip iStrip (rad)
  //
  

  Float_t delta =0.;
  if(iPlate==0)delta = -1. ;
  if(iPlate==1)delta = -0.5;
  if(iPlate==3)delta = +0.5;
  if(iPlate==4)delta = +1. ;

  Float_t r =(fgkRmin+fgkRmax)/2.+fHeights[iPlate][iStrip];
  Float_t z =(fgkRmin+2.)*TMath::Tan(fAngles[iPlate][iStrip]/kRaddeg)+delta;
  Float_t theta =0.5*TMath::Pi()-TMath::ATan(z/r);
  if(iPlate != 2){
  if(theta > 0.5*TMath::Pi() )theta+=fgkDprecCen;
  if(theta < 0.5*TMath::Pi() )theta-=fgkDprecCen;
  }
  return theta;

}
//_____________________________________________________________________________
void AliTOFGeometryV4::GetVolumePath(Int_t *ind, Char_t *path ) {
  //--------------------------------------------------------------------
  // This function returns the colume path of a given pad 
  //--------------------------------------------------------------------
  Int_t sector = ind[0];
  Char_t  string1[100];
  Char_t  string2[100];
  Char_t  string3[100];
  Char_t  string4[100];
  Int_t nstrB = NStripB();
  Int_t nstrC = NStripC();
  
  Int_t icopy=-1;
  
  if(sector<3){
    icopy=sector+1;
    sprintf(string1,"/ALIC_1/B077_1/B075_%i/BTO3_1",icopy);
  }
  else if(sector<11){
    //	icopy=sector-2;
    icopy=sector+3;
    sprintf(string1,"/ALIC_1/B077_1/B071_%i/BTO1_1",icopy);
  }
  else if(sector==11 || sector==12){
    icopy=sector-10;
    sprintf(string1,"/ALIC_1/B077_1/B074_%i/BTO2_1",icopy);
  }
  else {
    //	icopy=sector-4;
    icopy=sector-12;
    sprintf(string1,"/ALIC_1/B077_1/B071_%i/BTO1_1",icopy);
  }
  
  Int_t modnum=ind[1];
  Int_t istrip=ind[2];
  
  if( modnum ==0){
    sprintf(string2,"FTOC_1/FLTC_0");
    icopy= nstrC - istrip;
    sprintf(string3,"FSTR_%i",icopy);
  }    
  else if( modnum ==1){
    sprintf(string2,"FTOB_1/FLTB_0");
    icopy= nstrB - istrip;
      sprintf(string3,"FSTR_%i",icopy);
  }
  else if( modnum ==2){
    sprintf(string2,"FTOA_0/FLTA_0");
    icopy= istrip+1;
    sprintf(string3,"FSTR_%i",icopy);
  }
  else if( modnum ==3){
    sprintf(string2,"FTOB_2/FLTB_0");
    icopy= istrip+1;
    sprintf(string3,"FSTR_%i",icopy);
  }
  else if( modnum ==4){
    sprintf(string2,"FTOC_2/FLTC_0");
    icopy= istrip+1;
    sprintf(string3,"FSTR_%i",icopy);
  }


  Int_t padz = ind[3]+1; 
  Int_t padx = ind[4]+1;
  if(modnum==3 || modnum==4){
    padz = NpadZ() -ind[3];
    padx = NpadX() -ind[4];
  }
  sprintf(string4,"FSEN_0/FSEZ_%i/FSEX_%i",padz,padx);
  sprintf(path,"%s/%s/%s/%s",string1,string2,string3,string4); 

}

//_____________________________________________________________________________
void AliTOFGeometryV4::GetPos(Int_t *det, Float_t *pos) 
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
