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
#include "AliTRDrecPoint.h"

ClassImp(AliTRDgeometry)

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

}

//_____________________________________________________________________________
void AliTRDgeometry::Init()
{
  //
  // Initializes the geometry parameter
  //

  Int_t iplan;

  // The width of the chambers
  fCwidth[0]    =  99.6;
  fCwidth[1]    = 104.1;
  fCwidth[2]    = 108.5;
  fCwidth[3]    = 112.9;
  fCwidth[4]    = 117.4;
  fCwidth[5]    = 121.8;

  // The default pad dimensions
  fRowPadSize  = 4.5;
  fColPadSize  = 1.0;
  fTimeBinSize = 0.1;

  // The maximum number of pads
  // and the position of pad 0,0,0 
  // 
  // chambers seen from the top:
  //     +----------------------------+
  //     |                            |
  //     |                            |     ^
  //     |                            | rphi|
  //     |                            |     |
  //     |0                           |     | 
  //     +----------------------------+     +------>
  //                                             z 
  // chambers seen from the side:           ^
  //     +----------------------------+ time|
  //     |                            |     |
  //     |0                           |     |
  //     +----------------------------+     +------>
  //                                             z
  //                                             

  // The pad column (rphi-direction)  
  for (iplan = 0; iplan < kNplan; iplan++) {
    fColMax[iplan] = 1 + TMath::Nint((fCwidth[iplan] - 2. * kCcthick) 
                                                     / fColPadSize - 0.5);
    fCol0[iplan]   = -fCwidth[iplan]/2. + kCcthick;
  }

  // The time bucket
  fTimeMax = 1 + TMath::Nint(kDrThick / fTimeBinSize - 0.5);
  for (iplan = 0; iplan < kNplan; iplan++) {
    fTime0[iplan]  = kRmin + kCcframe/2. + kDrZpos - 0.5 * kDrThick
                           + iplan * (kCheight + kCspace);
  } 

}

//_____________________________________________________________________________
void AliTRDgeometry::CreateGeometry(Int_t *idtmed)
{
  //
  // Create the TRD geometry
  //
  // Author: Christoph Blume (C.Blume@gsi.de) 20/07/99
  //
  // The volumes:
  //    TRD1-3     (Air)   --- The TRD mother volumes for one sector. 
  //                           To be placed into the spaceframe.
  //
  //    UAFI(/M/O) (Al)    --- The aluminum frame of the inner(/middle/outer) chambers (readout)
  //    UCFI(/M/O) (C)     --- The carbon frame of the inner(/middle/outer) chambers 
  //                           (driftchamber + radiator)
  //    UAII(/M/O) (Air)   --- The inner part of the readout of the inner(/middle/outer) chambers
  //    UFII(/M/O) (Air)   --- The inner part of the chamner and radiator of the 
  //                           inner(/middle/outer) chambers
  //
  // The material layers in one chamber:
  //    UL01       (G10)   --- The gas seal of the radiator
  //    UL02       (CO2)   --- The gas in the radiator
  //    UL03       (PE)    --- The foil stack
  //    UL04       (Mylar) --- Entrance window to the driftvolume and HV-cathode
  //    UL05       (Xe)    --- The driftvolume
  //    UL06       (Xe)    --- The amplification region
  //    
  //    UL07       (Cu)    --- The pad plane
  //    UL08       (G10)   --- The Nomex honeycomb support structure
  //    UL09       (Cu)    --- FEE and signal lines
  //    UL10       (PE)    --- The cooling devices
  //    UL11       (Water) --- The cooling water

  const Int_t npar_cha = 3;

  Float_t par_dum[3];
  Float_t par_cha[npar_cha];

  Float_t xpos, ypos, zpos;

  // The aluminum frames - readout + electronics (Al)
  // The inner chambers
  gMC->Gsvolu("UAFI","BOX ",idtmed[1301-1],par_dum,0);
  // The middle chambers
  gMC->Gsvolu("UAFM","BOX ",idtmed[1301-1],par_dum,0);
  // The outer chambers
  gMC->Gsvolu("UAFO","BOX ",idtmed[1301-1],par_dum,0);

  // The inner part of the aluminum frames (Air)
  // The inner chambers
  gMC->Gsvolu("UAII","BOX ",idtmed[1302-1],par_dum,0);
  // The middle chambers
  gMC->Gsvolu("UAIM","BOX ",idtmed[1302-1],par_dum,0);
  // The outer chambers
  gMC->Gsvolu("UAIO","BOX ",idtmed[1302-1],par_dum,0);

  // The carbon frames - radiator + driftchamber (C)
  // The inner chambers
  gMC->Gsvolu("UCFI","BOX ",idtmed[1307-1],par_dum,0);
  // The middle chambers
  gMC->Gsvolu("UCFM","BOX ",idtmed[1307-1],par_dum,0);
  // The outer chambers
  gMC->Gsvolu("UCFO","BOX ",idtmed[1307-1],par_dum,0);

  // The inner part of the carbon frames (Air)
  // The inner chambers
  gMC->Gsvolu("UCII","BOX ",idtmed[1302-1],par_dum,0);
  // The middle chambers
  gMC->Gsvolu("UCIM","BOX ",idtmed[1302-1],par_dum,0);
  // The outer chambers
  gMC->Gsvolu("UCIO","BOX ",idtmed[1302-1],par_dum,0);

  // The material layers inside the chambers
  par_cha[0] = -1.;
  par_cha[1] = -1.;
  // G10 layer (radiator seal)
  par_cha[2] = kSeThick/2;
  gMC->Gsvolu("UL01","BOX ",idtmed[1313-1],par_cha,npar_cha);
  // CO2 layer (radiator)
  par_cha[2] = kRaThick/2;
  gMC->Gsvolu("UL02","BOX ",idtmed[1312-1],par_cha,npar_cha);
  // PE layer (radiator)
  par_cha[2] = kPeThick/2;
  gMC->Gsvolu("UL03","BOX ",idtmed[1303-1],par_cha,npar_cha);
  // Mylar layer (entrance window + HV cathode) 
  par_cha[2] = kMyThick/2;
  gMC->Gsvolu("UL04","BOX ",idtmed[1308-1],par_cha,npar_cha);
  // Xe/Isobutane layer (drift volume, sensitive) 
  par_cha[2] = kDrThick/2.;
  gMC->Gsvolu("UL05","BOX ",idtmed[1309-1],par_cha,npar_cha);
  // Xe/Isobutane layer (amplification volume, not sensitive)
  par_cha[2] = kAmThick/2.;
  gMC->Gsvolu("UL06","BOX ",idtmed[1309-1],par_cha,npar_cha);
  
  // Cu layer (pad plane)
  par_cha[2] = kCuThick/2;
  gMC->Gsvolu("UL07","BOX ",idtmed[1305-1],par_cha,npar_cha);
  // G10 layer (support structure)
  par_cha[2] = kSuThick/2;
  gMC->Gsvolu("UL08","BOX ",idtmed[1313-1],par_cha,npar_cha);
  // Cu layer (FEE + signal lines)
  par_cha[2] = kFeThick/2;
  gMC->Gsvolu("UL09","BOX ",idtmed[1305-1],par_cha,npar_cha);
  // PE layer (cooling devices)
  par_cha[2] = kCoThick/2;
  gMC->Gsvolu("UL10","BOX ",idtmed[1303-1],par_cha,npar_cha);
  // Water layer (cooling)
  par_cha[2] = kWaThick/2;
  gMC->Gsvolu("UL11","BOX ",idtmed[1314-1],par_cha,npar_cha);

  // Position the layers in the chambers
  xpos = 0;
  ypos = 0;

  // G10 layer (radiator seal)
  zpos = kSeZpos;
  gMC->Gspos("UL01",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL01",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL01",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // CO2 layer (radiator)
  zpos = kRaZpos;
  gMC->Gspos("UL02",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL02",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL02",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // PE layer (radiator)
  zpos = 0;
  gMC->Gspos("UL03",1,"UL02",xpos,ypos,zpos,0,"ONLY");
  // Mylar layer (entrance window + HV cathode)   
  zpos = kMyZpos;
  gMC->Gspos("UL04",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL04",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL04",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // Xe/Isobutane layer (drift volume) 
  zpos = kDrZpos;
  gMC->Gspos("UL05",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL05",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL05",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // Xe/Isobutane layer (amplification volume)
  zpos = kAmZpos;
  gMC->Gspos("UL06",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL06",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL06",3,"UCIO",xpos,ypos,zpos,0,"ONLY");

  // Cu layer (pad plane)
  zpos = kCuZpos;
  gMC->Gspos("UL07",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL07",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL07",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // G10 layer (support structure)
  zpos = kSuZpos;
  gMC->Gspos("UL08",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL08",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL08",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // Cu layer (FEE + signal lines)
  zpos = kFeZpos; 
  gMC->Gspos("UL09",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL09",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL09",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // PE layer (cooling devices)
  zpos = kCoZpos;
  gMC->Gspos("UL10",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL10",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL10",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // Water layer (cooling)
  zpos = kWaZpos;
  gMC->Gspos("UL11",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL11",1,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL11",1,"UAIO",xpos,ypos,zpos,0,"ONLY");

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::Local2Global(Int_t idet, Float_t *local, Float_t *global)
{
  //
  // Converts local pad-coordinates (row,col,time) into 
  // global ALICE reference frame coordinates (x,y,z)
  //

  Int_t        icham     = GetChamber(idet);    // Chamber info (0-4)
  Int_t        isect     = GetSector(idet);     // Sector info  (0-17)
  Int_t        iplan     = GetPlane(idet);      // Plane info   (0-5)

  return Local2Global(iplan,icham,isect,local,global);

}
 
//_____________________________________________________________________________
Bool_t AliTRDgeometry::Local2Global(Int_t iplan, Int_t icham, Int_t isect
                                  , Float_t *local, Float_t *global)
{
  //
  // Converts local pad-coordinates (row,col,time) into 
  // global ALICE reference frame coordinates (x,y,z)
  //

  Int_t        idet      = GetDetector(iplan,icham,isect); // Detector number

  Float_t      padRow    = local[0];                       // Pad Row position
  Float_t      padCol    = local[1];                       // Pad Column position
  Float_t      timeSlice = local[2];                       // Time "position"

  Float_t      row0      = GetRow0(iplan,icham,isect);
  Float_t      col0      = GetCol0(iplan);
  Float_t      time0     = GetTime0(iplan);

  Float_t      rot[3];

  // calculate (x,y,z) position in rotated chamber
  rot[0] = time0 + timeSlice * fTimeBinSize;
  rot[1] = col0  + padCol    * fColPadSize;
  rot[2] = row0  + padRow    * fRowPadSize;

  // Rotate back to original position
  return RotateBack(idet,rot,global);

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::Rotate(Int_t d, Float_t *pos, Float_t *rot)
{
  //
  // Rotates all chambers in the position of sector 0 and transforms
  // the coordinates in the ALICE restframe <pos> into the 
  // corresponding local frame <rot>.
  //

  Int_t   sector = GetSector(d);

  Float_t phi    = -2.0 * kPI /  (Float_t) kNsect * ((Float_t) sector + 0.5);

  rot[0] =  pos[0] * TMath::Cos(phi) + pos[1] * TMath::Sin(phi);
  rot[1] = -pos[0] * TMath::Sin(phi) + pos[1] * TMath::Cos(phi);
  rot[2] =  pos[2];

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::RotateBack(Int_t d, Float_t *rot, Float_t *pos)
{
  //
  // Rotates a chambers from the position of sector 0 into its
  // original position and transforms the corresponding local frame 
  // coordinates <rot> into the coordinates of the ALICE restframe <pos>.
  //

  Int_t   sector = GetSector(d);

  Float_t phi    =  2.0 * kPI /  (Float_t) kNsect * ((Float_t) sector + 0.5);

  pos[0] =  rot[0] * TMath::Cos(phi) + rot[1] * TMath::Sin(phi);
  pos[1] = -rot[0] * TMath::Sin(phi) + rot[1] * TMath::Cos(phi);
  pos[2] =  rot[2];

  return kTRUE;

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetDetector(Int_t p, Int_t c, Int_t s)
{
  //
  // Convert plane / chamber / sector into detector number
  //

  return (p + c * kNplan + s * kNplan * kNcham);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetPlane(Int_t d)
{
  //
  // Reconstruct the plane number from the detector number
  //

  return ((Int_t) (d % kNplan));

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetChamber(Int_t d)
{
  //
  // Reconstruct the chamber number from the detector number
  //

  return ((Int_t) (d % (kNplan * kNcham)) / kNplan);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetSector(Int_t d)
{
  //
  // Reconstruct the sector number from the detector number
  //

  return ((Int_t) (d / (kNplan * kNcham)));

}

//_____________________________________________________________________________
void AliTRDgeometry::GetGlobal(const AliRecPoint *p, TVector3 &pos, TMatrix &mat)
{
  // 
  // Returns the global coordinate and error matrix of a AliTRDrecPoint
  //

  GetGlobal(p,pos);

}

//_____________________________________________________________________________
void AliTRDgeometry::GetGlobal(const AliRecPoint *p, TVector3 &pos)
{
  // 
  // Returns the global coordinate and error matrix of a AliTRDrecPoint
  //

  Int_t detector = ((AliTRDrecPoint *) p)->GetDetector();

  Float_t global[3];
  Float_t local[3];
  local[0] = ((AliTRDrecPoint *) p)->GetLocalRow();
  local[1] = ((AliTRDrecPoint *) p)->GetLocalCol();
  local[2] = ((AliTRDrecPoint *) p)->GetLocalTime();

  if (Local2Global(detector,local,global)) {
    pos.SetX(global[0]);
    pos.SetY(global[1]);
    pos.SetZ(global[2]);
  }
  else {
    pos.SetX(0.0);
    pos.SetY(0.0);
    pos.SetZ(0.0);
  }

}
