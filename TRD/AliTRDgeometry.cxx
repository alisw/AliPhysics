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
Revision 1.22  2002/11/21 22:38:47  alibrary
Removing AliMC and AliMCProcess

Revision 1.21  2002/11/21 16:09:44  cblume
Change fgkRpadW to 1.0 cm for new pad plane

Revision 1.20  2002/10/31 17:45:35  cblume
New chamber geometry

Revision 1.19  2002/10/14 14:57:43  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.15.6.2  2002/07/24 10:09:30  alibrary
Updating VirtualMC

Revision 1.15.6.1  2002/06/10 15:28:58  hristov
Merged with v3-08-02

Revision 1.17  2002/04/05 13:20:12  cblume
Remove const for CreateGeometry

Revision 1.16  2002/03/28 14:59:07  cblume
Coding conventions

Revision 1.18  2002/06/12 09:54:35  cblume
Update of tracking code provided by Sergei

Revision 1.17  2002/04/05 13:20:12  cblume
Remove const for CreateGeometry

Revision 1.16  2002/03/28 14:59:07  cblume
Coding conventions

Revision 1.15  2002/02/11 14:21:16  cblume
Update of the geometry. Get rid of MANY

Revision 1.14  2001/11/06 17:19:41  cblume
Add detailed geometry and simple simulator

Revision 1.13  2001/08/02 08:30:45  cblume
Fix positions of cooling material

Revision 1.12  2001/05/21 16:45:47  hristov
Last minute changes (C.Blume)

Revision 1.11  2001/05/11 07:56:12  hristov
Consistent declarations needed on Alpha

Revision 1.10  2001/05/07 08:08:05  cblume
Update of TRD code

Revision 1.9  2001/03/27 12:48:33  cblume
Correct for volume overlaps

Revision 1.8  2001/03/13 09:30:35  cblume
Update of digitization. Moved digit branch definition to AliTRD

Revision 1.7  2001/02/14 18:22:26  cblume
Change in the geometry of the padplane

Revision 1.6  2000/11/01 14:53:20  cblume
Merge with TRD-develop

Revision 1.1.4.7  2000/10/16 01:16:53  cblume
Changed timebin 0 to be the one closest to the readout

Revision 1.1.4.6  2000/10/15 23:35:57  cblume
Include geometry constants as static member

Revision 1.1.4.5  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.4.4  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.1.4.3  2000/09/22 14:43:40  cblume
Allow the pad/timebin-dimensions to be changed after initialization

Revision 1.1.4.2  2000/09/18 13:37:01  cblume
Minor coding corrections

Revision 1.5  2000/10/02 21:28:19  fca
Removal of useless dependecies via forward declarations

Revision 1.4  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.3  2000/06/07 16:25:37  cblume
Try to remove compiler warnings on Sun and HP

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.1  2000/05/08 14:45:55  cblume
Bug fix in RotateBack(). Geometry update

Revision 1.4  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.3  2000/06/07 16:25:37  cblume
Try to remove compiler warnings on Sun and HP

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.1  2000/05/08 14:45:55  cblume
Bug fix in RotateBack(). Geometry update

Revision 1.1  2000/02/28 19:00:44  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD geometry class                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTRDgeometry.h"
#include "AliTRDparameter.h"

ClassImp(AliTRDgeometry)

//_____________________________________________________________________________

  //
  // The geometry constants
  //
  const Int_t   AliTRDgeometry::fgkNsect   = kNsect;
  const Int_t   AliTRDgeometry::fgkNplan   = kNplan;
  const Int_t   AliTRDgeometry::fgkNcham   = kNcham;
  const Int_t   AliTRDgeometry::fgkNdet    = kNdet;

  //
  // Dimensions of the detector
  //

  // Inner and outer radius of the mother volumes 
  const Float_t AliTRDgeometry::fgkRmin    = 294.0;
  const Float_t AliTRDgeometry::fgkRmax    = 368.0;

  // Upper and lower length of the mother volumes 
  const Float_t AliTRDgeometry::fgkZmax1   = 378.35; 
  const Float_t AliTRDgeometry::fgkZmax2   = 302.0; 

  // Parameter of the BTR mother volumes 
  const Float_t AliTRDgeometry::fgkSheight =  74.0; 
  const Float_t AliTRDgeometry::fgkSwidth1 =  99.613;
  const Float_t AliTRDgeometry::fgkSwidth2 = 125.707;
  const Float_t AliTRDgeometry::fgkSlenTR1 = 751.0;
  const Float_t AliTRDgeometry::fgkSlenTR2 = 313.5; 
  const Float_t AliTRDgeometry::fgkSlenTR3 = 159.5;  

  // The super module side plates
  const Float_t AliTRDgeometry::fgkSMpltT  = 0.2;
  const Float_t AliTRDgeometry::fgkSMgapT  = 0.5;  

  // Height of different chamber parts
  // Radiator
  const Float_t AliTRDgeometry::fgkCraH    =   4.8; 
  // Drift region
  const Float_t AliTRDgeometry::fgkCdrH    =   3.0;
  // Amplification region
  const Float_t AliTRDgeometry::fgkCamH    =   0.7;
  // Readout
  const Float_t AliTRDgeometry::fgkCroH    =   2.316;
  // Total height
  const Float_t AliTRDgeometry::fgkCH      = AliTRDgeometry::fgkCraH
                                           + AliTRDgeometry::fgkCdrH
                                           + AliTRDgeometry::fgkCamH
                                           + AliTRDgeometry::fgkCroH;  

  // Vertical spacing of the chambers
  const Float_t AliTRDgeometry::fgkVspace  =   1.784;

  // Horizontal spacing of the chambers
  const Float_t AliTRDgeometry::fgkHspace  =   2.0;

  // Thicknesses of different parts of the chamber frame
  // Lower aluminum frame
  const Float_t AliTRDgeometry::fgkCalT    =   0.3;
  // Lower G10 frame sides
  const Float_t AliTRDgeometry::fgkCclsT   =   0.3;
  // Lower G10 frame front
  const Float_t AliTRDgeometry::fgkCclfT   =   1.0;
  // Upper G10 frame
  const Float_t AliTRDgeometry::fgkCcuT    =   0.9;
  // Upper Al frame
  const Float_t AliTRDgeometry::fgkCauT    =   1.5;

  // Additional width of the readout chamber frames
  const Float_t AliTRDgeometry::fgkCroW    =   0.9;

  // Difference of outer chamber width and pad plane width
  //const Float_t AliTRDgeometry::fgkCpadW   =   1.0;
  const Float_t AliTRDgeometry::fgkCpadW   =   0.0;
  const Float_t AliTRDgeometry::fgkRpadW   =   1.0;

  //
  // Thickness of the the material layers
  //
  const Float_t AliTRDgeometry::fgkRaThick = 0.3646;  
  const Float_t AliTRDgeometry::fgkMyThick = 0.005;
  const Float_t AliTRDgeometry::fgkDrThick = AliTRDgeometry::fgkCdrH;    
  const Float_t AliTRDgeometry::fgkAmThick = AliTRDgeometry::fgkCamH;
  const Float_t AliTRDgeometry::fgkXeThick = AliTRDgeometry::fgkDrThick
                                           + AliTRDgeometry::fgkAmThick;
  const Float_t AliTRDgeometry::fgkCuThick = 0.001; 
  const Float_t AliTRDgeometry::fgkSuThick = 0.06; 
  const Float_t AliTRDgeometry::fgkFeThick = 0.0044; 
  const Float_t AliTRDgeometry::fgkCoThick = 0.02;
  const Float_t AliTRDgeometry::fgkWaThick = 0.02;

  //
  // Position of the material layers
  //
  const Float_t AliTRDgeometry::fgkRaZpos  = -1.50;
  const Float_t AliTRDgeometry::fgkMyZpos  =  0.895;
  const Float_t AliTRDgeometry::fgkDrZpos  =  2.4;
  const Float_t AliTRDgeometry::fgkAmZpos  =  0.0;
  const Float_t AliTRDgeometry::fgkCuZpos  = -0.9995;
  const Float_t AliTRDgeometry::fgkSuZpos  =  0.0000;
  const Float_t AliTRDgeometry::fgkFeZpos  =  0.0322;
  const Float_t AliTRDgeometry::fgkCoZpos  =  0.97;
  const Float_t AliTRDgeometry::fgkWaZpos  =  0.99;

//_____________________________________________________________________________
AliTRDgeometry::AliTRDgeometry():AliGeometry()
{
  //
  // AliTRDgeometry default constructor
  //

  Init();

}

//_____________________________________________________________________________
AliTRDgeometry::~AliTRDgeometry()
{
  //
  // AliTRDgeometry destructor
  //

}

//_____________________________________________________________________________
void AliTRDgeometry::Init()
{
  //
  // Initializes the geometry parameter
  //
  // The maximum number of pads
  // and the position of pad 0,0,0 
  // 
  // chambers seen from the top:
  //     +----------------------------+
  //     |                            |
  //     |                            |      ^
  //     |                            |  rphi|
  //     |                            |      |
  //     |0                           |      | 
  //     +----------------------------+      +------>
  //                                             z 
  // chambers seen from the side:            ^
  //     +----------------------------+ drift|
  //     |0                           |      |
  //     |                            |      |
  //     +----------------------------+      +------>
  //                                             z
  //                                             
  // IMPORTANT: time bin 0 is now the first one in the drift region 
  // closest to the readout !!!
  //

  Int_t icham;
  Int_t iplan;
  Int_t isect;

  // The outer width of the chambers
  fCwidth[0] =  94.8;
  fCwidth[1] =  99.3;
  fCwidth[2] = 103.7;
  fCwidth[3] = 108.1;
  fCwidth[4] = 112.6;
  fCwidth[5] = 117.0;

  // The outer lengths of the chambers
  // Includes the spacings between the chambers!
  Float_t length[kNplan][kNcham]   = { { 124.0, 124.0, 110.0, 124.0, 124.0 }
                                     , { 131.0, 131.0, 110.0, 131.0, 131.0 }
                                     , { 138.0, 138.0, 110.0, 138.0, 138.0 }
                                     , { 145.0, 145.0, 110.0, 145.0, 145.0 }
                                     , { 147.0, 147.0, 110.0, 147.0, 147.0 }
                                     , { 147.0, 147.0, 110.0, 147.0, 147.0 } };

  for (icham = 0; icham < kNcham; icham++) {
    for (iplan = 0; iplan < kNplan; iplan++) {
      fClength[iplan][icham]   = length[iplan][icham];
      fClengthPH[iplan][icham] = 0.0;
      fClengthRH[iplan][icham] = 0.0;
    }
  }

  // The rotation matrix elements
  Float_t phi = 0;
  for (isect = 0; isect < fgkNsect; isect++) {
    phi = -2.0 * TMath::Pi() /  (Float_t) fgkNsect * ((Float_t) isect + 0.5);
    fRotA11[isect] = TMath::Cos(phi);
    fRotA12[isect] = TMath::Sin(phi);
    fRotA21[isect] = TMath::Sin(phi);
    fRotA22[isect] = TMath::Cos(phi);
    phi = -1.0 * phi;
    fRotB11[isect] = TMath::Cos(phi);
    fRotB12[isect] = TMath::Sin(phi);
    fRotB21[isect] = TMath::Sin(phi);
    fRotB22[isect] = TMath::Cos(phi);
  }
 
}

//_____________________________________________________________________________
void AliTRDgeometry::CreateGeometry(Int_t* )
{
  //
  // Create TRD geometry
  //

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::Local2Global(Int_t idet, Float_t *local
                                   , Float_t *global
                                   , AliTRDparameter *par) const
{
  //
  // Converts local pad-coordinates (row,col,time) into 
  // global ALICE reference frame coordinates (x,y,z)
  //

  Int_t icham = GetChamber(idet);    // Chamber info (0-4)
  Int_t isect = GetSector(idet);     // Sector info  (0-17)
  Int_t iplan = GetPlane(idet);      // Plane info   (0-5)

  return Local2Global(iplan,icham,isect,local,global,par);

}
 
//_____________________________________________________________________________
Bool_t AliTRDgeometry::Local2Global(Int_t iplan, Int_t icham, Int_t isect
                                  , Float_t *local, Float_t *global
                                  , AliTRDparameter *par) const
{
  //
  // Converts local pad-coordinates (row,col,time) into 
  // global ALICE reference frame coordinates (x,y,z)
  //

  if (!par) {
    Error("Local2Global","No parameter defined\n");
    return kFALSE;
  }

  Int_t    idet      = GetDetector(iplan,icham,isect); // Detector number

  Float_t  padRow    = local[0]+0.5;                   // Pad Row position
  Float_t  padCol    = local[1]+0.5;                   // Pad Column position
  Float_t  timeSlice = local[2]+0.5;                   // Time "position"

  Float_t  row0      = par->GetRow0(iplan,icham,isect);
  Float_t  col0      = par->GetCol0(iplan);
  Float_t  time0     = par->GetTime0(iplan);

  Float_t  rot[3];

  // calculate (x,y,z) position in rotated chamber
  rot[0] = time0 - (timeSlice - par->GetTimeBefore()) 
         * par->GetTimeBinSize();
  rot[1] = col0  + padCol                    
         * par->GetColPadSize(iplan);
  rot[2] = row0  + padRow                    
         * par->GetRowPadSize(iplan,icham,isect);

  // Rotate back to original position
  return RotateBack(idet,rot,global);

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::Rotate(Int_t d, Float_t *pos, Float_t *rot) const
{
  //
  // Rotates all chambers in the position of sector 0 and transforms
  // the coordinates in the ALICE restframe <pos> into the 
  // corresponding local frame <rot>.
  //

  Int_t sector = GetSector(d);

  rot[0] =  pos[0] * fRotA11[sector] + pos[1] * fRotA12[sector];
  rot[1] = -pos[0] * fRotA21[sector] + pos[1] * fRotA22[sector];
  rot[2] =  pos[2];

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::RotateBack(Int_t d, Float_t *rot, Float_t *pos) const
{
  //
  // Rotates a chambers from the position of sector 0 into its
  // original position and transforms the corresponding local frame 
  // coordinates <rot> into the coordinates of the ALICE restframe <pos>.
  //

  Int_t sector = GetSector(d);

  pos[0] =  rot[0] * fRotB11[sector] + rot[1] * fRotB12[sector];
  pos[1] = -rot[0] * fRotB21[sector] + rot[1] * fRotB22[sector];
  pos[2] =  rot[2];

  return kTRUE;

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetDetectorSec(const Int_t p, const Int_t c) const
{
  //
  // Convert plane / chamber into detector number for one single sector
  //

  return (p + c * fgkNplan);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetDetector(const Int_t p, const Int_t c, const Int_t s) const
{
  //
  // Convert plane / chamber / sector into detector number
  //

  return (p + c * fgkNplan + s * fgkNplan * fgkNcham);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetPlane(const Int_t d) const
{
  //
  // Reconstruct the plane number from the detector number
  //

  return ((Int_t) (d % fgkNplan));

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetChamber(const Int_t d) const
{
  //
  // Reconstruct the chamber number from the detector number
  //

  return ((Int_t) (d % (fgkNplan * fgkNcham)) / fgkNplan);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetSector(const Int_t d) const
{
  //
  // Reconstruct the sector number from the detector number
  //

  return ((Int_t) (d / (fgkNplan * fgkNcham)));

}

//_____________________________________________________________________________
void AliTRDgeometry::SetOldGeometry()
{
  //
  // Use the old chamber lengths
  //

  Int_t icham;
  Int_t iplan;

  Float_t length[kNplan][kNcham]   = { { 123.5, 123.5, 110.0, 123.5, 123.5 }
				     , { 131.0, 131.0, 110.0, 131.0, 131.0 }
				     , { 134.5, 138.5, 110.0, 138.5, 134.5 }
				     , { 142.0, 146.0, 110.0, 146.0, 142.0 }
				     , { 142.0, 153.0, 110.0, 153.0, 142.0 }
                                     , { 134.0, 160.5, 110.0, 160.5, 134.0 } };

  for (icham = 0; icham < kNcham; icham++) {
    for (iplan = 0; iplan < kNplan; iplan++) {
      fClength[iplan][icham]   = length[iplan][icham];
    }
  }

}
