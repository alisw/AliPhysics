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

const Float_t AliTOFGeometry::fgkXPad     = 2.5;    // Pad size in the x direction (cm)
const Float_t AliTOFGeometry::fgkZPad     = 3.5;    // Pad size in the z direction (cm)

const Float_t AliTOFGeometry::fgkStripLength = 122.;// Strip Length (rho X phi direction) (cm)

const Float_t AliTOFGeometry::fgkSigmaForTail1= 2.; //Sig1 for simulation of TDC tails 
const Float_t AliTOFGeometry::fgkSigmaForTail2= 0.5;//Sig2 for simulation of TDC tails

const Float_t AliTOFGeometry::fgkTdcBin = 24.4;     // time-window for the TDC bins [ps]

//_____________________________________________________________________________
AliTOFGeometry::AliTOFGeometry()
{
  //
  // AliTOFGeometry default constructor
  //

  fNStripC     = 19;  // number of strips in C type module 
  fZlenA    = 106.0;  // length (cm) of the A module
  fZlenB    = 141.0;  // length (cm) of the B module
  fZlenC    = 177.5;  // length (cm) of the C module
  fMaxhZtof = 371.5;  // Max half z-size of TOF (cm)

  fxTOF     = 371.; // Inner radius of the TOF for Reconstruction (cm)
  fRmin     = 370.; // Inner radius of the TOF (cm)
  fRmax     = 399.; // Outer radius of the TOF (cm)

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
 
  fPhiSec   = 360./kNSectors;

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

   for (Int_t iplate = 0; iplate < kNPlates; iplate++) {
     for (Int_t istrip = 0; istrip < kMaxNstrip; istrip++) {
       fAngles[iplate][istrip]   = kangles[iplate][istrip];
       fHeights[iplate][istrip]  = kheights[iplate][istrip];
     }
   }

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
