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

const Float_t AliTOFGeometry::fgkSigmaForTail1= 2.; //Sig1 for simulation of TDC tails 
const Float_t AliTOFGeometry::fgkSigmaForTail2= 0.5;//Sig2 for simulation of TDC tails

const Float_t AliTOFGeometry::fgkTdcBin = 24.4;     // time-window for the TDC bins [ps]

//_____________________________________________________________________________
AliTOFGeometry::AliTOFGeometry()
{
  //
  // AliTOFGeometry default constructor
  //

  kNStripC     = 20;  // number of strips in C type module 
  kMaxNstrip   = 20;  // Max. number of strips 
  kZlenA    = 106.0;  // length (cm) of the A module
  kZlenB    = 141.0;  // length (cm) of the B module
  kZlenC    = 177.5;  // length (cm) of the C module
  kMaxhZtof = 371.5;  // Max half z-size of TOF (cm)
  kStripLength = 122.;// Strip Length (rho X phi direction) (cm)

  fgkxTOF     = 371.; // Inner radius of the TOF for Reconstruction (cm)
  fgkRmin     = 370.; // Inner radius of the TOF (cm)
  fgkRmax     = 399.; // Outer radius of the TOF (cm)

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
