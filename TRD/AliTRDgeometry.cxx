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

#include "AliMC.h"

#include "AliTRDgeometry.h"
#include "AliTRDrecPoint.h"
#include "AliMC.h"

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
  const Float_t AliTRDgeometry::fgkRmin    = 294.0;
  const Float_t AliTRDgeometry::fgkRmax    = 368.0;

  const Float_t AliTRDgeometry::fgkZmax1   = 378.35; 
  const Float_t AliTRDgeometry::fgkZmax2   = 302.0; 

  const Float_t AliTRDgeometry::fgkSheight =  74.0; 
  const Float_t AliTRDgeometry::fgkSwidth1 =  99.613;
  const Float_t AliTRDgeometry::fgkSwidth2 = 125.707;
  const Float_t AliTRDgeometry::fgkSlenTR1 = 751.0;
  const Float_t AliTRDgeometry::fgkSlenTR2 = 313.5; 
  const Float_t AliTRDgeometry::fgkSlenTR3 = 159.5;  

  const Float_t AliTRDgeometry::fgkCheight =  11.0;  
  const Float_t AliTRDgeometry::fgkCspace  =   1.6;
  const Float_t AliTRDgeometry::fgkCathick =   1.0; 
  const Float_t AliTRDgeometry::fgkCcthick =   1.0;
  const Float_t AliTRDgeometry::fgkCaframe =   2.675; 
  const Float_t AliTRDgeometry::fgkCcframe = AliTRDgeometry::fgkCheight 
                                           - AliTRDgeometry::fgkCaframe;

  //
  // Thickness of the the material layers
  //
  const Float_t AliTRDgeometry::fgkSeThick = 0.02;  
  const Float_t AliTRDgeometry::fgkRaThick = 4.8;  
  const Float_t AliTRDgeometry::fgkPeThick = 0.20;    
  const Float_t AliTRDgeometry::fgkMyThick = 0.005;
  const Float_t AliTRDgeometry::fgkXeThick = 3.5;
  const Float_t AliTRDgeometry::fgkDrThick = 3.0;
  const Float_t AliTRDgeometry::fgkAmThick = AliTRDgeometry::fgkXeThick 
                                           - AliTRDgeometry::fgkDrThick;
  const Float_t AliTRDgeometry::fgkCuThick = 0.001; 
  const Float_t AliTRDgeometry::fgkSuThick = 0.06; 
  const Float_t AliTRDgeometry::fgkFeThick = 0.0044; 
  const Float_t AliTRDgeometry::fgkCoThick = 0.02;
  const Float_t AliTRDgeometry::fgkWaThick = 0.01;

  //
  // Position of the material layers
  //
  const Float_t AliTRDgeometry::fgkSeZpos  = -4.1525; 
  const Float_t AliTRDgeometry::fgkRaZpos  = -1.7425;
  const Float_t AliTRDgeometry::fgkPeZpos  =  0.0000;
  const Float_t AliTRDgeometry::fgkMyZpos  =  0.6600;
  const Float_t AliTRDgeometry::fgkDrZpos  =  2.1625;
  const Float_t AliTRDgeometry::fgkAmZpos  =  4.1125;
  const Float_t AliTRDgeometry::fgkCuZpos  = -1.3370; 
  const Float_t AliTRDgeometry::fgkSuZpos  =  0.0000;
  const Float_t AliTRDgeometry::fgkFeZpos  =  1.3053;
  const Float_t AliTRDgeometry::fgkCoZpos  =  1.3175;
  const Float_t AliTRDgeometry::fgkWaZpos  =  1.3325; 

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

  Int_t isect;

  // The width of the chambers
  fCwidth[0] =  99.6;
  fCwidth[1] = 104.1;
  fCwidth[2] = 108.5;
  fCwidth[3] = 112.9;
  fCwidth[4] = 117.4;
  fCwidth[5] = 121.8;

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
  // IMPORTANT: time bin 0 is now the one closest to the readout !!!
  //

  // The pad column (rphi-direction)  
  SetColPadSize(1.0);

  // The time bucket
  SetTimeBinSize(0.1);

  // The rotation matrix elements
  Float_t phi = 0;
  for (isect = 0; isect < fgkNsect; isect++) {
    phi = -2.0 * kPI /  (Float_t) fgkNsect * ((Float_t) isect + 0.5);
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
void AliTRDgeometry::SetColPadSize(Float_t size)
{
  //
  // Redefines the pad size in column direction
  //

  fColPadSize = size;
  for (Int_t iplan = 0; iplan < fgkNplan; iplan++) {
    fColMax[iplan] = 1 + TMath::Nint((fCwidth[iplan] - 2. * fgkCcthick) 
                                                     / fColPadSize - 0.5);
    fCol0[iplan]   = -fCwidth[iplan]/2. + fgkCcthick;
  }

}

//_____________________________________________________________________________
void AliTRDgeometry::SetTimeBinSize(Float_t size)
{
  //
  // Redefines the time bin size
  //

  fTimeBinSize = size;
  fTimeMax     = 1 + TMath::Nint(fgkDrThick / fTimeBinSize - 0.5);
  for (Int_t iplan = 0; iplan < fgkNplan; iplan++) {
    fTime0[iplan]  = fgkRmin + fgkCcframe/2. + fgkDrZpos + 0.5 * fgkDrThick
                             + iplan * (fgkCheight + fgkCspace);
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

  const Int_t kNparCha = 3;

  Float_t parDum[3];
  Float_t parCha[kNparCha];

  Float_t xpos, ypos, zpos;

  // The aluminum frames - readout + electronics (Al)
  // The inner chambers
  gMC->Gsvolu("UAFI","BOX ",idtmed[1301-1],parDum,0);
  // The middle chambers
  gMC->Gsvolu("UAFM","BOX ",idtmed[1301-1],parDum,0);
  // The outer chambers
  gMC->Gsvolu("UAFO","BOX ",idtmed[1301-1],parDum,0);

  // The inner part of the aluminum frames (Air)
  // The inner chambers
  gMC->Gsvolu("UAII","BOX ",idtmed[1302-1],parDum,0);
  // The middle chambers
  gMC->Gsvolu("UAIM","BOX ",idtmed[1302-1],parDum,0);
  // The outer chambers
  gMC->Gsvolu("UAIO","BOX ",idtmed[1302-1],parDum,0);

  // The carbon frames - radiator + driftchamber (C)
  // The inner chambers
  gMC->Gsvolu("UCFI","BOX ",idtmed[1307-1],parDum,0);
  // The middle chambers
  gMC->Gsvolu("UCFM","BOX ",idtmed[1307-1],parDum,0);
  // The outer chambers
  gMC->Gsvolu("UCFO","BOX ",idtmed[1307-1],parDum,0);

  // The inner part of the carbon frames (Air)
  // The inner chambers
  gMC->Gsvolu("UCII","BOX ",idtmed[1302-1],parDum,0);
  // The middle chambers
  gMC->Gsvolu("UCIM","BOX ",idtmed[1302-1],parDum,0);
  // The outer chambers
  gMC->Gsvolu("UCIO","BOX ",idtmed[1302-1],parDum,0);

  // The material layers inside the chambers
  parCha[0] = -1.;
  parCha[1] = -1.;
  // G10 layer (radiator seal)
  parCha[2] = fgkSeThick/2;
  gMC->Gsvolu("UL01","BOX ",idtmed[1313-1],parCha,kNparCha);
  // CO2 layer (radiator)
  parCha[2] = fgkRaThick/2;
  gMC->Gsvolu("UL02","BOX ",idtmed[1312-1],parCha,kNparCha);
  // PE layer (radiator)
  parCha[2] = fgkPeThick/2;
  gMC->Gsvolu("UL03","BOX ",idtmed[1303-1],parCha,kNparCha);
  // Mylar layer (entrance window + HV cathode) 
  parCha[2] = fgkMyThick/2;
  gMC->Gsvolu("UL04","BOX ",idtmed[1308-1],parCha,kNparCha);
  // Xe/Isobutane layer (drift volume, sensitive) 
  parCha[2] = fgkDrThick/2.;
  gMC->Gsvolu("UL05","BOX ",idtmed[1309-1],parCha,kNparCha);
  // Xe/Isobutane layer (amplification volume, not sensitive)
  parCha[2] = fgkAmThick/2.;
  gMC->Gsvolu("UL06","BOX ",idtmed[1309-1],parCha,kNparCha);
  
  // Cu layer (pad plane)
  parCha[2] = fgkCuThick/2;
  gMC->Gsvolu("UL07","BOX ",idtmed[1305-1],parCha,kNparCha);
  // G10 layer (support structure)
  parCha[2] = fgkSuThick/2;
  gMC->Gsvolu("UL08","BOX ",idtmed[1313-1],parCha,kNparCha);
  // Cu layer (FEE + signal lines)
  parCha[2] = fgkFeThick/2;
  gMC->Gsvolu("UL09","BOX ",idtmed[1305-1],parCha,kNparCha);
  // PE layer (cooling devices)
  parCha[2] = fgkCoThick/2;
  gMC->Gsvolu("UL10","BOX ",idtmed[1303-1],parCha,kNparCha);
  // Water layer (cooling)
  parCha[2] = fgkWaThick/2;
  gMC->Gsvolu("UL11","BOX ",idtmed[1314-1],parCha,kNparCha);

  // Position the layers in the chambers
  xpos = 0;
  ypos = 0;

  // G10 layer (radiator seal)
  zpos = fgkSeZpos;
  gMC->Gspos("UL01",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL01",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL01",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // CO2 layer (radiator)
  zpos = fgkRaZpos;
  gMC->Gspos("UL02",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL02",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL02",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // PE layer (radiator)
  zpos = 0;
  gMC->Gspos("UL03",1,"UL02",xpos,ypos,zpos,0,"ONLY");
  // Mylar layer (entrance window + HV cathode)   
  zpos = fgkMyZpos;
  gMC->Gspos("UL04",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL04",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL04",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // Xe/Isobutane layer (drift volume) 
  zpos = fgkDrZpos;
  gMC->Gspos("UL05",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL05",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL05",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // Xe/Isobutane layer (amplification volume)
  zpos = fgkAmZpos;
  gMC->Gspos("UL06",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL06",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL06",3,"UCIO",xpos,ypos,zpos,0,"ONLY");

  // Cu layer (pad plane)
  zpos = fgkCuZpos;
  gMC->Gspos("UL07",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL07",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL07",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // G10 layer (support structure)
  zpos = fgkSuZpos;
  gMC->Gspos("UL08",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL08",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL08",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // Cu layer (FEE + signal lines)
  zpos = fgkFeZpos; 
  gMC->Gspos("UL09",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL09",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL09",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // PE layer (cooling devices)
  zpos = fgkCoZpos;
  gMC->Gspos("UL10",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL10",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL10",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // Water layer (cooling)
  zpos = fgkWaZpos;
  gMC->Gspos("UL11",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL11",1,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL11",1,"UAIO",xpos,ypos,zpos,0,"ONLY");

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::Local2Global(Int_t idet, Float_t *local, Float_t *global) const
{
  //
  // Converts local pad-coordinates (row,col,time) into 
  // global ALICE reference frame coordinates (x,y,z)
  //

  Int_t icham = GetChamber(idet);    // Chamber info (0-4)
  Int_t isect = GetSector(idet);     // Sector info  (0-17)
  Int_t iplan = GetPlane(idet);      // Plane info   (0-5)

  return Local2Global(iplan,icham,isect,local,global);

}
 
//_____________________________________________________________________________
Bool_t AliTRDgeometry::Local2Global(Int_t iplan, Int_t icham, Int_t isect
                                  , Float_t *local, Float_t *global) const
{
  //
  // Converts local pad-coordinates (row,col,time) into 
  // global ALICE reference frame coordinates (x,y,z)
  //

  Int_t    idet      = GetDetector(iplan,icham,isect); // Detector number

  Float_t  padRow    = local[0];                       // Pad Row position
  Float_t  padCol    = local[1];                       // Pad Column position
  Float_t  timeSlice = local[2];                       // Time "position"

  Float_t  row0      = GetRow0(iplan,icham,isect);
  Float_t  col0      = GetCol0(iplan);
  Float_t  time0     = GetTime0(iplan);

  Float_t  rot[3];

  // calculate (x,y,z) position in rotated chamber
  rot[0] = time0 + timeSlice * fTimeBinSize;
  rot[1] = col0  + padCol    * fColPadSize;
  rot[2] = row0  + padRow    * fRowPadSize;

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
Int_t AliTRDgeometry::GetDetector(Int_t p, Int_t c, Int_t s) const
{
  //
  // Convert plane / chamber / sector into detector number
  //

  return (p + c * fgkNplan + s * fgkNplan * fgkNcham);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetPlane(Int_t d) const
{
  //
  // Reconstruct the plane number from the detector number
  //

  return ((Int_t) (d % fgkNplan));

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetChamber(Int_t d) const
{
  //
  // Reconstruct the chamber number from the detector number
  //

  return ((Int_t) (d % (fgkNplan * fgkNcham)) / fgkNplan);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetSector(Int_t d) const
{
  //
  // Reconstruct the sector number from the detector number
  //

  return ((Int_t) (d / (fgkNplan * fgkNcham)));

}

//_____________________________________________________________________________
void AliTRDgeometry::GetGlobal(const AliRecPoint *p, TVector3 &pos
                             , TMatrix &mat) const
{
  // 
  // Returns the global coordinate and error matrix of a AliTRDrecPoint
  //

  GetGlobal(p,pos);
  mat.Zero();

}

//_____________________________________________________________________________
void AliTRDgeometry::GetGlobal(const AliRecPoint *p, TVector3 &pos) const
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
