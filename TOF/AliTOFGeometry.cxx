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

#include <stdlib.h>
#include <Riostream.h>
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TOF Geometry class                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliConst.h"
#include "AliTOFGeometry.h"

ClassImp(AliTOFGeometry)

const Int_t AliTOFGeometry::fgkTimeDiff   = 25000;  // Min signal separation (ps)

const Float_t AliTOFGeometry::fgkxTOF     = 371.;   // Inner radius of the TOF for Reconstruction (cm)
const Float_t AliTOFGeometry::fgkRmin     = 370.;   // Inner radius of the TOF (cm)
const Float_t AliTOFGeometry::fgkRmax     = 399;    // Outer radius of the TOF (cm)
const Float_t AliTOFGeometry::fgkZlenA    = 106.0;  // length (cm) of the A module
const Float_t AliTOFGeometry::fgkZlenB    = 141.0;  // length (cm) of the B module
const Float_t AliTOFGeometry::fgkZlenC    = 177.5;  // length (cm) of the C module
const Float_t AliTOFGeometry::fgkXPad     = 2.5;    // Pad size in the x direction (cm)
const Float_t AliTOFGeometry::fgkZPad     = 3.5;    // Pad size in the z direction (cm)
const Float_t AliTOFGeometry::fgkMaxhZtof = 371.5;  // Max half z-size of TOF (cm)
const Float_t AliTOFGeometry::fgkStripLength = 122.;// Strip Length (rho X phi direction) (cm)
const Float_t AliTOFGeometry::fgkDeadBndX = 1.0;    // Dead Boundaries of a Strip along X direction (length) (cm)
const Float_t AliTOFGeometry::fgkDeadBndZ = 1.5;    // Dead Boundaries of a Strip along Z direction (width) (cm)
const Float_t AliTOFGeometry::fgkOverSpc = 15.3;    // Space available for sensitive layers in radial direction (cm)


const Float_t AliTOFGeometry::fgkSigmaForTail1= 2.;//Sig1 for simulation of TDC tails 
const Float_t AliTOFGeometry::fgkSigmaForTail2= 0.5;//Sig2 for simulation of TDC tails
const Float_t AliTOFGeometry::fgkSpeedOfLight = 0.299792458;// c (10^9 m/s)
const Float_t AliTOFGeometry::fgkPionMass     = 0.13957;// pion mass (Gev/c^2)
const Float_t AliTOFGeometry::fgkKaonMass     = 0.49368;// kaon mass (Gev/c^2)
const Float_t AliTOFGeometry::fgkProtonMass   = 0.93827;// proton mass (Gev/c^2)
const Float_t AliTOFGeometry::fgkElectronMass = 0.00051;// electron mass (Gev/c^2)
const Float_t AliTOFGeometry::fgkMuonMass     = 0.10566;// muon mass (Gev/c^2)


const Float_t AliTOFGeometry::fgkDprecMin = 0.0000075;//num.prec.tolerance on Thmin 
const Float_t AliTOFGeometry::fgkDprecMax = 0.0000100;//num.prec.tolerance on Thma 
const Float_t AliTOFGeometry::fgkDprecCen = 0.0000005;//num.prec.tolerance on <Theta> 

const Float_t AliTOFGeometry::fgkTdcBin = 24.4; // time-window for the TDC bins [ps]

//_____________________________________________________________________________
AliTOFGeometry::AliTOFGeometry()
{
  //
  // AliTOFGeometry default constructor
  //
  Init();

}

//_____________________________________________________________________________
AliTOFGeometry::~AliTOFGeometry()
{
  //
  // AliTOFGeometry destructor
  //

}
//_____________________________________________________________________________
void AliTOFGeometry::Init()
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
      fAngles[iplate][istrip]   = kangles[iplate][istrip];
      fHeights[iplate][istrip]  = kheights[iplate][istrip];
    }
  }

  fPhiSec   = 360./kNSectors;
}




//_____________________________________________________________________________
Float_t AliTOFGeometry::DistanceToPad(Int_t *det, Float_t *pos, Float_t *dist3d) 
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
  Float_t angle   = phi*kRaddeg-( Int_t (kRaddeg*phi/20.) + 0.5)*fPhiSec;
  Float_t xs = radius*TMath::Cos(angle/kRaddeg);
  Float_t ys = radius*TMath::Sin(angle/kRaddeg);
  Float_t zs = z;

  // Do the same for the selected pad

  Float_t g[3];
  GetPos(det,g);

  Float_t padRadius = TMath::Sqrt(g[0]*g[0]+g[1]*g[1]);
  Float_t padPhi=TMath::ATan2(g[1],g[0]);	
  if(padPhi<0) padPhi=2.*TMath::Pi()+padPhi;
  //  Get the local angle in the sector philoc
  Float_t padAngle   = padPhi*kRaddeg-( Int_t (padPhi*kRaddeg/20.)+ 0.5) * fPhiSec; 
  Float_t padxs = padRadius*TMath::Cos(padAngle/kRaddeg);
  Float_t padys = padRadius*TMath::Sin(padAngle/kRaddeg);
  Float_t padzs = g[2];
  
  //Now move to local pad coordinate frame. Translate:
  
  Float_t xt = xs-padxs;
  Float_t yt = ys-padys;
  Float_t zt = zs-padzs;
  //Now Rotate:
  
  Float_t alpha = GetAngles(det[1],det[2]);
  Float_t xr = xt*TMath::Cos(alpha/kRaddeg)+zt*TMath::Sin(alpha/kRaddeg);
  Float_t yr = yt;
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
Bool_t AliTOFGeometry::IsInsideThePad(Int_t *det, Float_t *pos) 
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
  Float_t angle   = phi*kRaddeg-( Int_t (kRaddeg*phi/20.) + 0.5) *fPhiSec;
  Float_t xs = radius*TMath::Cos(angle/kRaddeg);
  Float_t ys = radius*TMath::Sin(angle/kRaddeg);
  Float_t zs = z;

  // Do the same for the selected pad

  Float_t g[3];
  GetPos(det,g);

  Float_t padRadius = TMath::Sqrt(g[0]*g[0]+g[1]*g[1]);
  Float_t padPhi=TMath::ATan2(g[1],g[0]);	
  if(padPhi<0) padPhi=2.*TMath::Pi()+padPhi;
  //  Get the local angle in the sector philoc
  Float_t padAngle   = padPhi*kRaddeg-( Int_t (padPhi*kRaddeg/20.)+ 0.5) * fPhiSec; 
  Float_t padxs = padRadius*TMath::Cos(padAngle/kRaddeg);
  Float_t padys = padRadius*TMath::Sin(padAngle/kRaddeg);
  Float_t padzs = g[2];
  
  //Now move to local pad coordinate frame. Translate:
  
  Float_t xt = xs-padxs;
  Float_t yt = ys-padys;
  Float_t zt = zs-padzs;
  //Now Rotate:
  
  Float_t alpha = GetAngles(det[1],det[2]);
  Float_t xr = xt*TMath::Cos(alpha/kRaddeg)+zt*TMath::Sin(alpha/kRaddeg);
  Float_t yr = yt;
  Float_t zr = -xt*TMath::Sin(alpha/kRaddeg)+zt*TMath::Cos(alpha/kRaddeg);

  if(TMath::Abs(xr)<=0.75 && TMath::Abs(yr)<= (fgkXPad*0.5) && TMath::Abs(zr)<= (fgkZPad*0.5))
    isInside=true; 
  return isInside;

}

//_____________________________________________________________________________
void AliTOFGeometry::GetPos(Int_t *det, Float_t *pos) 
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
void AliTOFGeometry::GetDetID( Float_t *pos, Int_t *det) 
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
Float_t AliTOFGeometry::GetX(Int_t *det) 
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
Float_t AliTOFGeometry::GetY(Int_t *det) 
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
Float_t AliTOFGeometry::GetZ(Int_t *det) 
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
Int_t AliTOFGeometry::GetSector(Float_t *pos) 
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
Int_t AliTOFGeometry::GetPadX(Float_t *pos) 
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
Int_t AliTOFGeometry::GetPlate(Float_t *pos) 
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
Int_t AliTOFGeometry::GetStrip(Float_t *pos) 
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
Int_t AliTOFGeometry::GetPadZ(Float_t *pos) 
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
Float_t AliTOFGeometry::GetMinPlateTheta(Int_t iPlate) 
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
Float_t AliTOFGeometry::GetMaxPlateTheta(Int_t iPlate) 
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
Float_t  AliTOFGeometry::GetMaxStripTheta(Int_t iPlate, Int_t iStrip) 
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
Float_t  AliTOFGeometry::GetMinStripTheta(Int_t iPlate, Int_t iStrip) 
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
Float_t  AliTOFGeometry::GetStripTheta(Int_t iPlate, Int_t iStrip) 
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



