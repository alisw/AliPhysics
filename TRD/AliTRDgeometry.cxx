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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD geometry class                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TError.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TGeoMatrix.h>


#include "AliRunLoader.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"

#include "AliAlignObj.h"
#include "AliAlignObjAngles.h"

#include "AliRun.h"
#include "AliTRD.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"

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
  const Float_t AliTRDgeometry::fgkSMpltT  =   0.2;
  const Float_t AliTRDgeometry::fgkSMgapT  =   0.5;  

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
  
  const Double_t AliTRDgeometry::fgkTime0Base = Rmin() + CraHght() + CdrHght() + CamHght()/2.;
  const Float_t  AliTRDgeometry::fgkTime0[6]  = { fgkTime0Base + 0 * (Cheight() + Cspace()), 
                                                  fgkTime0Base + 1 * (Cheight() + Cspace()), 
                                                  fgkTime0Base + 2 * (Cheight() + Cspace()), 
                                                  fgkTime0Base + 3 * (Cheight() + Cspace()), 
                                                  fgkTime0Base + 4 * (Cheight() + Cspace()), 
                                                  fgkTime0Base + 5 * (Cheight() + Cspace()) };

//_____________________________________________________________________________
AliTRDgeometry::AliTRDgeometry():AliGeometry()
{
  //
  // AliTRDgeometry default constructor
  //
  fMatrixArray =0;
  fMatrixCorrectionArray= 0;
  Init();
}

//_____________________________________________________________________________
AliTRDgeometry::~AliTRDgeometry()
{
  //
  // AliTRDgeometry destructor
  //
  delete fMatrixArray;
  delete fMatrixCorrectionArray;
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
  //
  // Changed with the introduction of 
  // the new layer 0. The old layer 6
  // is removed.
  fCwidth[0] =  90.4;
  fCwidth[1] =  94.8;
  fCwidth[2] =  99.3;
  fCwidth[3] = 103.7;
  fCwidth[4] = 108.1;
  fCwidth[5] = 112.6;
  // Old layer 6
  // fCwidth[5] = 117.0;

  // The outer lengths of the chambers
  // Includes the spacings between the chambers!
  // Changed with the introduction of 
  // the new layer 0. The old layer 6
  // is removed.
  Float_t length[kNplan][kNcham]   = { { 124.0, 124.0, 110.0, 124.0, 124.0 }
				     , { 124.0, 124.0, 110.0, 124.0, 124.0 }
                                     , { 131.0, 131.0, 110.0, 131.0, 131.0 }
                                     , { 138.0, 138.0, 110.0, 138.0, 138.0 }
                                     , { 145.0, 145.0, 110.0, 145.0, 145.0 }
				     , { 147.0, 147.0, 110.0, 147.0, 147.0 } };
  // Old layer 6
  //                                 , { 147.0, 147.0, 110.0, 147.0, 147.0 } };

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
Bool_t AliTRDgeometry::Local2Global(Int_t idet, Double_t *local
                                   , Double_t *global) const
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
                                  , Double_t *local, Double_t *global) const
{
  //
  // Converts local pad-coordinates (row,col,time) into 
  // global ALICE reference frame coordinates (x,y,z)
  //

  AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();
  if (!commonParam)
    return kFALSE;

  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration)
    return kFALSE;  
  
  AliTRDpadPlane *padPlane = commonParam->GetPadPlane(iplan,icham);

  // calculate (x,y,z) position in rotated chamber
  Int_t    row       = ((Int_t) local[0]);
  Int_t    col       = ((Int_t) local[1]);
  Float_t  timeSlice = local[2] + 0.5;
  Float_t  time0     = GetTime0(iplan);

  Int_t idet = GetDetector(iplan, icham, isect);

  Double_t  rot[3];
  rot[0] = time0 - (timeSlice - calibration->GetT0(idet, col, row))
      * calibration->GetVdrift(idet, col, row)/calibration->GetSamplingFrequency();
  rot[1] = padPlane->GetColPos(col) - 0.5 * padPlane->GetColSize(col);
  rot[2] = padPlane->GetRowPos(row) - 0.5 * padPlane->GetRowSize(row);

  // Rotate back to original position
  return RotateBack(idet,rot,global);

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::Global2Local(Int_t mode, Double_t *local, Double_t *global
                                   , Int_t* index) const
{
  //
  // Converts local pad-coordinates (row,col,time) into 
  // global ALICE reference frame coordinates (x,y,z)
  //
  // index[0] = plane number
  // index[1] = chamber number
  // index[2] = sector number
  //
  // mode=0  - local coordinate in y, z,             x - rotated global   
  // mode=2  - local coordinate in pad, and pad row, x - rotated global
  //

  //Int_t    idet    = GetDetector(iplan,icham,isect); // Detector number
  Int_t    idet      = GetDetector(index[0],index[1],index[2]); // Detector number
  RotateBack(idet,global,local);
  if (mode==0) return kTRUE;
  //
  //  Float_t  row0      = par->GetRow0(iplan,icham,isect);
  //Float_t  col0      = par->GetCol0(iplan);
  //Float_t  time0     = GetTime0(iplan);
  //
  // mode 1 to be implemented later
  // calculate (x,y,z) position in time bin pad row pad
  //
  //rot[0] = time0 - (timeSlice - par->GetTimeBefore()) 
  //       * par->GetDriftVelocity()/par->GetSamplingFrequency();
  //rot[1] = col0  + padCol                    
  //       * par->GetColPadSize(iplan);
  //rot[2] = row0  + padRow                    
  //       * par->GetRowPadSize(iplan,icham,isect);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::Global2Detector(Double_t global[3], Int_t index[3])
{
  //  
  //  Find detector for given global point - Ideal geometry 
  //  
  //
  // input    = global position
  // output   = index
  // index[0] = plane number
  // index[1] = chamber number
  // index[2] = sector number
  //
  //
  // Find sector
  //
  Float_t fi = TMath::ATan2(global[1],global[0]);
  if (fi<0) fi += 2*TMath::Pi();
  index[2] = fgkNsect-1-TMath::Nint((fi - GetAlpha()/2.)/GetAlpha());
  //
  // Find plane
  //
  Float_t locx = global[0] * fRotA11[index[2]] + global[1] * fRotA12[index[2]];  
  index[0] = 0;
  Float_t max = locx - GetTime0(0);
  for (Int_t iplane=1; iplane<fgkNplan;iplane++){
    Float_t dist = TMath::Abs(locx - GetTime0(iplane));
    if (dist < max){
      index[0] = iplane;
      max = dist;
    }
  }
  //
  // Find chamber
  //
  if (TMath::Abs(global[2]) < 0.5*GetChamberLength(index[0],2)){
    index[1]=2;
  }else{
    Double_t localZ = global[2];
    if (global[2]>0){
      localZ -= 0.5*(GetChamberLength(index[0],2)+GetChamberLength(index[0],1));
      index[1] = (TMath::Abs(localZ) < 0.5*GetChamberLength(index[0],3)) ? 1:0;
    }
    else{
      localZ += 0.5*(GetChamberLength(index[0],2)+GetChamberLength(index[0],3));
      index[1] = (TMath::Abs(localZ) < 0.5*GetChamberLength(index[0],1)) ? 3:4;
    }
  }  
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliTRDgeometry::Rotate(Int_t d, Double_t *pos, Double_t *rot) const
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
Bool_t AliTRDgeometry::RotateBack(Int_t d, Double_t *rot, Double_t *pos) const
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
Int_t AliTRDgeometry::GetDetectorSec(Int_t p, Int_t c)
{
  //
  // Convert plane / chamber into detector number for one single sector
  //

  return (p + c * fgkNplan);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetDetector(Int_t p, Int_t c, Int_t s)
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
AliTRDgeometry* AliTRDgeometry::GetGeometry(AliRunLoader* runLoader)
{
  //
  // load the geometry from the galice file
  //

  if (!runLoader) runLoader = AliRunLoader::GetRunLoader();
  if (!runLoader) {
    ::Error("AliTRDgeometry::GetGeometry", "No run loader");
    return NULL;
  }

  TDirectory* saveDir = gDirectory;
  runLoader->CdGAFile();

  // Try from the galice.root file
  AliTRDgeometry* geom = (AliTRDgeometry*) gDirectory->Get("TRDgeometry");

  if (!geom) {
    // It is not in the file, try to get it from gAlice, 
    // which corresponds to the run loader 
    AliTRD * trd = (AliTRD*)runLoader->GetAliRun()->GetDetector("TRD");
    geom = trd->GetGeometry();
  }
  if (!geom) ::Error("AliTRDgeometry::GetGeometry", "Geometry not found");

  saveDir->cd();
  return geom;
}


//_____________________________________________________________________________
Bool_t   AliTRDgeometry::ReadGeoMatrices(){
  //
  // Read geo matrices from current gGeoManager for each TRD sector
  //

  //
  // fMatrixArray - 
  //

  //
  // fMatrixCorrectionArray - 
  //


  if (!gGeoManager) return kFALSE;
  fMatrixArray = new TObjArray(kNdet); 
  fMatrixCorrectionArray = new TObjArray(kNdet);
  fMatrixGeo   = new TObjArray(kNdet);
  AliAlignObjAngles o;
  //
  for (Int_t iLayer = AliAlignObj::kTRD1; iLayer <= AliAlignObj::kTRD6; iLayer++) {
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(iLayer); iModule++) {
      UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iModule);
      const char *path = AliAlignObj::GetVolPath(volid);
      if (!gGeoManager->cd(path)) return kFALSE;      
      TGeoHMatrix* m = gGeoManager->GetCurrentMatrix();
      Int_t     iLayerTRD    = iLayer-AliAlignObj::kTRD1;
      Int_t     isector      = Nsect()-1-(iModule/Ncham());
      Int_t     ichamber     = Ncham()-1-(iModule%Ncham());
      Int_t     lid          = GetDetector(iLayerTRD,ichamber,isector);    
      //
      //
      //
      // local geo system z-x-y  to x-y--z 
      //
      fMatrixGeo->AddAt(new TGeoHMatrix(*m),lid);
      
      TGeoRotation mchange; 
      mchange.RotateY(90); mchange.RotateX(90);
      //
      TGeoHMatrix gMatrix(mchange.Inverse());
      gMatrix.MultiplyLeft(m);
      fMatrixArray->AddAt(new TGeoHMatrix(gMatrix),lid); 
      //
      //  Cluster transformation matrix
      //
      TGeoHMatrix  rotMatrix(mchange.Inverse());
      rotMatrix.MultiplyLeft(m);
      Double_t sectorAngle = 20.*(isector%18)+10;
      TGeoHMatrix  rotSector;
      rotSector.RotateZ(sectorAngle);
      rotMatrix.MultiplyLeft(&rotSector);      
      //
      fMatrixCorrectionArray->AddAt(new TGeoHMatrix(rotMatrix),lid);       
    }    
  }
  return kTRUE;
}


