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


#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TGeoMatrix.h>

#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliAlignObj.h"
#include "AliAlignObjAngles.h"
#include "AliRun.h"

#include "AliTRD.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"

ClassImp(AliTRDgeometry)

//_____________________________________________________________________________

  //
  // The geometry constants
  //
  const Int_t    AliTRDgeometry::fgkNsect     = kNsect;
  const Int_t    AliTRDgeometry::fgkNplan     = kNplan;
  const Int_t    AliTRDgeometry::fgkNcham     = kNcham;
  const Int_t    AliTRDgeometry::fgkNdet      = kNdet;

  //
  // Dimensions of the detector
  //

  // Parameter of the BTRD mother volumes 
  const Float_t  AliTRDgeometry::fgkSheight   =  77.9; 
  const Float_t  AliTRDgeometry::fgkSwidth1   =  94.881; 
  const Float_t  AliTRDgeometry::fgkSwidth2   = 122.353;
  const Float_t  AliTRDgeometry::fgkSlength   = 751.0;

  // The super module side plates
  const Float_t  AliTRDgeometry::fgkSMpltT    =   0.2;

  // Height of different chamber parts
  // Radiator
  const Float_t  AliTRDgeometry::fgkCraH      =   4.8; 
  // Drift region
  const Float_t  AliTRDgeometry::fgkCdrH      =   3.0;
  // Amplification region
  const Float_t  AliTRDgeometry::fgkCamH      =   0.7;
  // Readout
  const Float_t  AliTRDgeometry::fgkCroH      =   2.316;
  // Total height
  const Float_t  AliTRDgeometry::fgkCH        = AliTRDgeometry::fgkCraH
                                              + AliTRDgeometry::fgkCdrH
                                              + AliTRDgeometry::fgkCamH
                                              + AliTRDgeometry::fgkCroH;  

  // Vertical spacing of the chambers
  const Float_t  AliTRDgeometry::fgkVspace    =   1.784;
  // Horizontal spacing of the chambers
  const Float_t  AliTRDgeometry::fgkHspace    =   2.0;
  // Radial distance of the first ROC to the outer plates of the SM
  const Float_t  AliTRDgeometry::fgkVrocsm    =   1.2;

  // Thicknesses of different parts of the chamber frame
  // Lower aluminum frame
  const Float_t  AliTRDgeometry::fgkCalT      =   0.4;
  // Lower Wacosit frame sides
  const Float_t  AliTRDgeometry::fgkCclsT     =   0.21;
  // Lower Wacosit frame front
  const Float_t  AliTRDgeometry::fgkCclfT     =   1.0;
  // Thickness of glue around radiator
  const Float_t  AliTRDgeometry::fgkCglT      =   0.25;
  // Upper Wacosit frame
  const Float_t  AliTRDgeometry::fgkCcuT      =   0.9;
  // Al frame of back panel
  const Float_t  AliTRDgeometry::fgkCauT      =   1.5;
  // Additional Al of the lower chamber frame
  const Float_t  AliTRDgeometry::fgkCalW      =   1.11;

  // Additional width of the readout chamber frames
  const Float_t  AliTRDgeometry::fgkCroW      =   0.9;

  // Difference of outer chamber width and pad plane width
  const Float_t  AliTRDgeometry::fgkCpadW     =   0.0;
  const Float_t  AliTRDgeometry::fgkRpadW     =   1.0;

  //
  // Thickness of the the material layers
  //
  const Float_t  AliTRDgeometry::fgkMyThick   = 0.005;
  const Float_t  AliTRDgeometry::fgkRaThick   = 0.3233;  
  const Float_t  AliTRDgeometry::fgkDrThick   = AliTRDgeometry::fgkCdrH;    
  const Float_t  AliTRDgeometry::fgkAmThick   = AliTRDgeometry::fgkCamH;
  const Float_t  AliTRDgeometry::fgkXeThick   = AliTRDgeometry::fgkDrThick
                                              + AliTRDgeometry::fgkAmThick;
  const Float_t  AliTRDgeometry::fgkWrThick   = 0.0002;
  const Float_t  AliTRDgeometry::fgkCuThick   = 0.0072; 
  const Float_t  AliTRDgeometry::fgkGlThick   = 0.05;
  const Float_t  AliTRDgeometry::fgkSuThick   = 0.0919; 
  const Float_t  AliTRDgeometry::fgkRcThick   = 0.0058;
  const Float_t  AliTRDgeometry::fgkRpThick   = 0.0632;
  const Float_t  AliTRDgeometry::fgkRoThick   = 0.0028;

  //
  // Position of the material layers
  //
  const Float_t  AliTRDgeometry::fgkRaZpos    =  0.0;
  const Float_t  AliTRDgeometry::fgkDrZpos    =  2.4;
  const Float_t  AliTRDgeometry::fgkAmZpos    =  0.0;
  const Float_t  AliTRDgeometry::fgkWrZpos    =  0.0;
  const Float_t  AliTRDgeometry::fgkCuZpos    = -0.9995;
  const Float_t  AliTRDgeometry::fgkGlZpos    = -0.5; 
  const Float_t  AliTRDgeometry::fgkSuZpos    =  0.0;
  const Float_t  AliTRDgeometry::fgkRcZpos    =  1.04;
  const Float_t  AliTRDgeometry::fgkRpZpos    =  1.0;
  const Float_t  AliTRDgeometry::fgkRoZpos    =  1.05;

  const Int_t    AliTRDgeometry::fgkMCMmax    = 16;   
  const Int_t    AliTRDgeometry::fgkMCMrow    = 4;   
  const Int_t    AliTRDgeometry::fgkROBmaxC0  = 6; 
  const Int_t    AliTRDgeometry::fgkROBmaxC1  = 8; 
  const Int_t    AliTRDgeometry::fgkADCmax    = 21;   
  const Int_t    AliTRDgeometry::fgkTBmax     = 60;   
  const Int_t    AliTRDgeometry::fgkPadmax    = 18;   
  const Int_t    AliTRDgeometry::fgkColmax    = 144;
  const Int_t    AliTRDgeometry::fgkRowmaxC0  = 12;
  const Int_t    AliTRDgeometry::fgkRowmaxC1  = 16;

  const Double_t AliTRDgeometry::fgkTime0Base = 300.65;
  const Float_t  AliTRDgeometry::fgkTime0[6]  = { fgkTime0Base + 0 * (Cheight() + Cspace()) 
                                                , fgkTime0Base + 1 * (Cheight() + Cspace()) 
                                                , fgkTime0Base + 2 * (Cheight() + Cspace()) 
                                                , fgkTime0Base + 3 * (Cheight() + Cspace()) 
                                                , fgkTime0Base + 4 * (Cheight() + Cspace()) 
                                                , fgkTime0Base + 5 * (Cheight() + Cspace())};

//_____________________________________________________________________________
AliTRDgeometry::AliTRDgeometry()
  :AliGeometry()
  ,fMatrixArray(0)
  ,fMatrixCorrectionArray(0)
  ,fMatrixGeo(0)

{
  //
  // AliTRDgeometry default constructor
  //

  Init();

}

//_____________________________________________________________________________
AliTRDgeometry::AliTRDgeometry(const AliTRDgeometry &g)
  :AliGeometry(g)
  ,fMatrixArray(g.fMatrixArray)
  ,fMatrixCorrectionArray(g.fMatrixCorrectionArray)
  ,fMatrixGeo(g.fMatrixGeo)
{
  //
  // AliTRDgeometry copy constructor
  //

  Init();

}

//_____________________________________________________________________________
AliTRDgeometry::~AliTRDgeometry()
{
  //
  // AliTRDgeometry destructor
  //

  if (fMatrixArray) {
    delete fMatrixArray;
    fMatrixArray           = 0;
  }

  if (fMatrixCorrectionArray) {
    delete fMatrixCorrectionArray;
    fMatrixCorrectionArray = 0;
  }

}

//_____________________________________________________________________________
AliTRDgeometry &AliTRDgeometry::operator=(const AliTRDgeometry &g)
{
  //
  // Assignment operator
  //

  if (this != &g) {
    Init();
  }

  return *this;

}

//_____________________________________________________________________________
void AliTRDgeometry::Init()
{
  //
  // Initializes the geometry parameter
  //

  Int_t icham;
  Int_t iplan;
  Int_t isect;

  // The outer width of the chambers
  fCwidth[0] =  90.4;
  fCwidth[1] =  94.8;
  fCwidth[2] =  99.3;
  fCwidth[3] = 103.7;
  fCwidth[4] = 108.1;
  fCwidth[5] = 112.6;

  // The outer lengths of the chambers
  // Includes the spacings between the chambers!
  Float_t length[kNplan][kNcham]   = { { 124.0, 124.0, 110.0, 124.0, 124.0 }
				     , { 124.0, 124.0, 110.0, 124.0, 124.0 }
                                     , { 131.0, 131.0, 110.0, 131.0, 131.0 }
                                     , { 138.0, 138.0, 110.0, 138.0, 138.0 }
                                     , { 145.0, 145.0, 110.0, 145.0, 145.0 }
				     , { 147.0, 147.0, 110.0, 147.0, 147.0 } };

  for (icham = 0; icham < kNcham; icham++) {
    for (iplan = 0; iplan < kNplan; iplan++) {
      fClength[iplan][icham] = length[iplan][icham];
    }
  }

  // The rotation matrix elements
  Float_t phi = 0.0;
  for (isect = 0; isect < fgkNsect; isect++) {
    phi = 2.0 * TMath::Pi() /  (Float_t) fgkNsect * ((Float_t) isect + 0.5);
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

  for (isect = 0; isect < fgkNsect; isect++) {
    SetSMstatus(isect,1);
  }
 
}

//_____________________________________________________________________________
void AliTRDgeometry::CreateGeometry(Int_t *idtmed)
{
  //
  // Create the TRD geometry without hole
  //
  //
  // Names of the TRD volumina (xx = detector number):
  //
  //      Volume (Air) wrapping the readout chamber components
  //        UTxx    includes: UAxx, UDxx, UFxx, UUxx
  //
  //      Volume (Air) wrapping the services (fee + cooling)
  //        UUxx    the services volume has been reduced by 7.42 mm
  //                in order to allow shifts in radial direction
  //
  //      Lower part of the readout chambers (drift volume + radiator)
  //
  //        UAxx    Aluminum frames                 (Al)
  //        UBxx    Wacosit frames                  (C)
  //        UXxx    Glue around radiator            (Epoxy)
  //        UCxx    Inner volumes                   (Air)
  //        UZxx    Additional aluminum ledges      (Al)
  //
  //      Upper part of the readout chambers (readout plane + fee)
  //
  //        UDxx    Wacosit frames of amp. region   (C)
  //        UExx    Inner volumes of the frame      (Air)
  //        UFxx    Aluminum frame of back panel    (Al)
  //        UGxx    Inner volumes of the back panel (Air)
  //
  //      Inner material layers
  //
  //        UHxx    Radiator                        (Rohacell)
  //        UJxx    Drift volume                    (Xe/CO2)
  //        UKxx    Amplification volume            (Xe/CO2)
  //        UWxx    Wire plane                      (Cu)
  //        ULxx    Pad plane                       (Cu)
  //        UYxx    Glue layer                      (Epoxy)
  //        UMxx    Support structure               (Rohacell)
  //        UNxx    ROB base material               (C)
  //        UOxx    ROB copper                      (Cu)
  //        UVxx    ROB other materials             (Cu)
  //

  const Int_t kNparTrd = 4;
  const Int_t kNparCha = 3;

  Float_t xpos;
  Float_t ypos;
  Float_t zpos;

  Float_t parTrd[kNparTrd];
  Float_t parCha[kNparCha];

  Char_t  cTagV[6];
  Char_t  cTagM[5];

  // The TRD mother volume for one sector (Air), full length in z-direction
  // Provides material for side plates of super module
  parTrd[0] = fgkSwidth1/2.0;
  parTrd[1] = fgkSwidth2/2.0;
  parTrd[2] = fgkSlength/2.0;
  parTrd[3] = fgkSheight/2.0;
  gMC->Gsvolu("UTR1","TRD1",idtmed[1302-1],parTrd,kNparTrd);

  // The outer aluminum plates of the super module (Al)
  parTrd[0] = fgkSwidth1/2.0;
  parTrd[1] = fgkSwidth2/2.0;
  parTrd[2] = fgkSlength/2.0;
  parTrd[3] = fgkSheight/2.0;
  gMC->Gsvolu("UTS1","TRD1",idtmed[1301-1],parTrd,kNparTrd);

  // The inner part of the TRD mother volume for one sector (Air), 
  // full length in z-direction
  parTrd[0] = fgkSwidth1/2.0 - fgkSMpltT;
  parTrd[1] = fgkSwidth2/2.0 - fgkSMpltT;
  parTrd[2] = fgkSlength/2.0;
  parTrd[3] = fgkSheight/2.0 - fgkSMpltT;
  gMC->Gsvolu("UTI1","TRD1",idtmed[1302-1],parTrd,kNparTrd);

  for (Int_t icham = 0; icham < kNcham; icham++) {
    for (Int_t iplan = 0; iplan < kNplan; iplan++) {  

      Int_t iDet = GetDetectorSec(iplan,icham);

      // The lower part of the readout chambers (drift volume + radiator) 
      // The aluminum frames 
      sprintf(cTagV,"UA%02d",iDet);
      parCha[0] = fCwidth[iplan]/2.0;
      parCha[1] = fClength[iplan][icham]/2.0 - fgkHspace/2.0;
      parCha[2] = fgkCraH/2.0 + fgkCdrH/2.0;
      fChamberUAboxd[iDet][0] = parCha[0];
      fChamberUAboxd[iDet][1] = parCha[1];
      fChamberUAboxd[iDet][2] = parCha[2];
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
      // The additional aluminum on the frames
      // This part has not the correct postion but is just supposed to
      // represent the missing material. The correct from of the L-shaped
      // profile would not fit into the alignable volume. 
      sprintf(cTagV,"UZ%02d",iDet);
      parCha[0] = fgkCroW/2.0;
      parCha[1] = fClength[iplan][icham]/2.0 - fgkHspace/2.0;
      parCha[2] = fgkCalW/2.0;
      fChamberUAboxd[iDet][0] = fChamberUAboxd[iDet][0] + fgkCroW;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
      // The Wacosit frames 
      sprintf(cTagV,"UB%02d",iDet);
      parCha[0] = fCwidth[iplan]/2.0 - fgkCalT; 
      parCha[1] = -1.0;
      parCha[2] = -1.0;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
      // The glue around the radiator
      sprintf(cTagV,"UX%02d",iDet);
      parCha[0] = fCwidth[iplan]/2.0 - fgkCalT - fgkCclsT; 
      parCha[1] = fClength[iplan][icham]/2.0 - fgkHspace/2.0 - fgkCclfT;
      parCha[2] = fgkCraH/2.0;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1311-1],parCha,kNparCha);
      // The inner part of radiator (air)
      sprintf(cTagV,"UC%02d",iDet);
      parCha[0] = fCwidth[iplan]/2.0 - fgkCalT - fgkCclsT - fgkCglT; 
      parCha[1] = fClength[iplan][icham]/2.0 - fgkHspace/2.0 - fgkCclfT - fgkCglT;
      parCha[2] = -1.0;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);

      // The upper part of the readout chambers (amplification volume)
      // The Wacosit frames
      sprintf(cTagV,"UD%02d",iDet);
      parCha[0] = fCwidth[iplan]/2.0 + fgkCroW;
      parCha[1] = fClength[iplan][icham]/2.0 - fgkHspace/2.0;
      parCha[2] = fgkCamH/2.0;
      fChamberUDboxd[iDet][0] = parCha[0];
      fChamberUDboxd[iDet][1] = parCha[1];
      fChamberUDboxd[iDet][2] = parCha[2];
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
      // The inner part of the Wacosit frame (air)
      sprintf(cTagV,"UE%02d",iDet);
      parCha[0] = fCwidth[iplan]/2.0 + fgkCroW - fgkCcuT; 
      parCha[1] = fClength[iplan][icham]/2.0 - fgkHspace/2.0 - fgkCcuT;
      parCha[2] = -1.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);

      // The support structure (pad plane, back panel, readout boards)
      // The aluminum frames
      sprintf(cTagV,"UF%02d",iDet);
      parCha[0] = fCwidth[iplan]/2.0 + fgkCroW;
      parCha[1] = fClength[iplan][icham]/2.0 - fgkHspace/2.0;
      parCha[2] = fgkCroH/2.0;
      fChamberUFboxd[iDet][0] = parCha[0];
      fChamberUFboxd[iDet][1] = parCha[1];
      fChamberUFboxd[iDet][2] = parCha[2];
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
      // The inner part of the aluminum frames
      sprintf(cTagV,"UG%02d",iDet);
      parCha[0] = fCwidth[iplan]/2.0 + fgkCroW - fgkCauT; 
      parCha[1] = fClength[iplan][icham]/2.0 - fgkHspace/2.0 - fgkCauT;
      parCha[2] = -1.0;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);

      // The material layers inside the chambers
      // Rohacell layer (radiator)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkRaThick/2.0;
      sprintf(cTagV,"UH%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1315-1],parCha,kNparCha);
      // Xe/Isobutane layer (drift volume) 
      parCha[0] = fCwidth[iplan]/2.0 - fgkCalT - fgkCclsT;
      parCha[1] = fClength[iplan][icham]/2.0 - fgkHspace/2.0 - fgkCclfT;
      parCha[2] = fgkDrThick/2.0;
      sprintf(cTagV,"UJ%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1309-1],parCha,kNparCha);
      // Xe/Isobutane layer (amplification volume)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkAmThick/2.0;
      sprintf(cTagV,"UK%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1309-1],parCha,kNparCha);  
      // Cu layer (wire plane)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkWrThick/2.0;
      sprintf(cTagV,"UW%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1303-1],parCha,kNparCha);
      // Cu layer (pad plane)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkCuThick/2.0;
      sprintf(cTagV,"UL%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1305-1],parCha,kNparCha);
      // Epoxy layer (glue)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkGlThick/2.0;
      sprintf(cTagV,"UY%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1311-1],parCha,kNparCha);
      // G10 layer (support structure / honeycomb)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkSuThick/2.0;
      sprintf(cTagV,"UM%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1310-1],parCha,kNparCha);
      // G10 layer (PCB readout board)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkRpThick/2;
      sprintf(cTagV,"UN%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1313-1],parCha,kNparCha);
      // Cu layer (traces in readout board)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkRcThick/2.0;
      sprintf(cTagV,"UO%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1306-1],parCha,kNparCha);
      // Cu layer (other material on in readout board)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkRoThick/2.0;
      sprintf(cTagV,"UV%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1304-1],parCha,kNparCha);

      // Position the layers in the chambers
      xpos = 0.0;
      ypos = 0.0;
      // Lower part
      // Rohacell layer (radiator)
      zpos = fgkRaZpos;
      sprintf(cTagV,"UH%02d",iDet);
      sprintf(cTagM,"UC%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Xe/Isobutane layer (drift volume) 
      zpos = fgkDrZpos;
      sprintf(cTagV,"UJ%02d",iDet);
      sprintf(cTagM,"UB%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Upper part
      // Xe/Isobutane layer (amplification volume)
      zpos = fgkAmZpos;
      sprintf(cTagV,"UK%02d",iDet);
      sprintf(cTagM,"UE%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Cu layer (wire plane inside amplification volume)
      zpos = fgkWrZpos; 
      sprintf(cTagV,"UW%02d",iDet);
      sprintf(cTagM,"UK%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Readout part + support plane
      // Cu layer (pad plane)
      zpos = fgkCuZpos; 
      sprintf(cTagV,"UL%02d",iDet);
      sprintf(cTagM,"UG%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Epoxy layer (glue)
      zpos = fgkGlZpos; 
      sprintf(cTagV,"UY%02d",iDet);
      sprintf(cTagM,"UG%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // G10 layer (support structure)
      zpos = fgkSuZpos;
      sprintf(cTagV,"UM%02d",iDet);
      sprintf(cTagM,"UG%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // G10 layer (PCB readout board)
      zpos = fgkRpZpos;
      sprintf(cTagV,"UN%02d",iDet);
      sprintf(cTagM,"UG%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Cu layer (traces in readout board)
      zpos = fgkRcZpos;
      sprintf(cTagV,"UO%02d",iDet);
      sprintf(cTagM,"UG%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Cu layer (other materials on readout board)
      zpos = fgkRoZpos;
      sprintf(cTagV,"UV%02d",iDet);
      sprintf(cTagM,"UG%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");

      // Position the inner volumes of the chambers in the frames
      xpos = 0.0;
      ypos = 0.0;
      // The inner part of the radiator
      zpos = 0.0;
      sprintf(cTagV,"UC%02d",iDet);
      sprintf(cTagM,"UX%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // The glue around the radiator
      zpos = fgkCraH/2.0 - fgkCdrH/2.0 - fgkCraH/2.0;
      sprintf(cTagV,"UX%02d",iDet);
      sprintf(cTagM,"UB%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // The lower Wacosit frame inside the aluminum frame
      zpos = 0.0;
      sprintf(cTagV,"UB%02d",iDet);
      sprintf(cTagM,"UA%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // The inside of the upper Wacosit frame
      zpos = 0.0;
      sprintf(cTagV,"UE%02d",iDet);
      sprintf(cTagM,"UD%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // The inside of the upper aluminum frame
      zpos = 0.0;
      sprintf(cTagV,"UG%02d",iDet);
      sprintf(cTagM,"UF%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");      

      // Position the frames of the chambers in the TRD mother volume
      xpos  = 0.0;
      ypos  = - fClength[iplan][0] - fClength[iplan][1] - fClength[iplan][2]/2.0;
      for (Int_t ic = 0; ic < icham; ic++) {
        ypos += fClength[iplan][ic];        
      }
      ypos += fClength[iplan][icham]/2.0;
      zpos  = fgkVrocsm + fgkSMpltT + fgkCraH/2.0 + fgkCdrH/2.0 - fgkSheight/2.0 
            + iplan * (fgkCH + fgkVspace);
      // The lower aluminum frame, radiator + drift region
      sprintf(cTagV,"UA%02d",iDet);      
      fChamberUAorig[iDet][0] = xpos;
      fChamberUAorig[iDet][1] = ypos;
      fChamberUAorig[iDet][2] = zpos;
      // The upper G10 frame, amplification region
      sprintf(cTagV,"UD%02d",iDet);
      zpos += fgkCamH/2.0 + fgkCraH/2.0 + fgkCdrH/2.0;      
      fChamberUDorig[iDet][0] = xpos;
      fChamberUDorig[iDet][1] = ypos;
      fChamberUDorig[iDet][2] = zpos;
      // The upper aluminum frame
      sprintf(cTagV,"UF%02d",iDet);
      zpos += fgkCroH/2.0 + fgkCamH/2.0;      
      fChamberUForig[iDet][0] = xpos;
      fChamberUForig[iDet][1] = ypos;
      fChamberUForig[iDet][2] = zpos;

    }
  }

  // Create the volumes of the super module frame
  CreateFrame(idtmed);

  // Create the volumes of the services
  CreateServices(idtmed);
  
  for (Int_t icham = 0; icham < kNcham; icham++) {
    for (Int_t iplan = 0; iplan < kNplan; iplan++) {  
      GroupChamber(iplan,icham,idtmed);
    }
  }
  
  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  gMC->Gspos("UTI1",1,"UTS1",xpos,ypos,zpos,0,"ONLY");

  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  gMC->Gspos("UTS1",1,"UTR1",xpos,ypos,zpos,0,"ONLY");

  // Put the TRD volumes into the space frame mother volumes
  // if enabled via status flag
  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  for (Int_t isect = 0; isect < kNsect; isect++) {
    if (fSMstatus[isect]) {
      sprintf(cTagV,"BTRD%d",isect);
      gMC->Gspos("UTR1",1,cTagV,xpos,ypos,zpos,0,"ONLY");
    }
  }

}

//_____________________________________________________________________________
void AliTRDgeometry::CreateFrame(Int_t *idtmed)
{
  //
  // Create the geometry of the frame of the supermodule
  //
  // Names of the TRD services volumina
  //
  //        USRL    Support rails for the chambers (Al)
  //        USxx    Support cross bars between the chambers (Al)
  //        USHx    Horizontal connection between the cross bars (Al)
  //        USLx    Long corner ledges (Al)
  //

  Int_t   iplan = 0;

  Float_t xpos  = 0.0;
  Float_t ypos  = 0.0;
  Float_t zpos  = 0.0;

  Char_t  cTagV[5];
  Char_t  cTagM[5];

  // The rotation matrices
  const Int_t kNmatrix = 4;
  Int_t   matrix[kNmatrix];
  gMC->Matrix(matrix[0], 100.0,   0.0,  90.0,  90.0,  10.0,   0.0);
  gMC->Matrix(matrix[1],  80.0,   0.0,  90.0,  90.0,  10.0, 180.0);
  gMC->Matrix(matrix[2],  90.0,   0.0,   0.0,   0.0,  90.0,  90.0);
  gMC->Matrix(matrix[3],  90.0, 180.0,   0.0, 180.0,  90.0,  90.0);

  //
  // The chamber support rails
  //

  const Float_t kSRLwid  = 2.00;
  const Float_t kSRLhgt  = 2.3;
  const Float_t kSRLdst  = 1.0;
  const Int_t   kNparSRL = 3;
  Float_t parSRL[kNparSRL];
  parSRL[0] = kSRLwid   /2.0;
  parSRL[1] = fgkSlength/2.0;
  parSRL[2] = kSRLhgt   /2.0;
  gMC->Gsvolu("USRL","BOX ",idtmed[1301-1],parSRL,kNparSRL);

  xpos  = 0.0;
  ypos  = 0.0;
  zpos  = 0.0;
  for (iplan = 0; iplan < kNplan; iplan++) {
    xpos  = fCwidth[iplan]/2.0 + kSRLwid/2.0 + kSRLdst;
    ypos  = 0.0;
    zpos  = fgkVrocsm + fgkSMpltT + fgkCraH + fgkCdrH + fgkCamH 
          - fgkSheight/2.0  
          + iplan * (fgkCH + fgkVspace);
    gMC->Gspos("USRL",iplan+1         ,"UTI1", xpos,ypos,zpos,0,"ONLY");
    gMC->Gspos("USRL",iplan+1+  kNplan,"UTI1",-xpos,ypos,zpos,0,"ONLY");
  }

  //
  // The cross bars between the chambers
  //

  const Float_t kSCBwid  = 1.0;
  const Float_t kSCBthk  = 2.0;
  const Float_t kSCHhgt  = 0.3;

  const Int_t   kNparSCB = 3;
  Float_t parSCB[kNparSCB];
  parSCB[1] = kSCBwid/2.0;
  parSCB[2] = fgkCH  /2.0 + fgkVspace/2.0 - kSCHhgt;

  const Int_t   kNparSCI = 3;
  Float_t parSCI[kNparSCI];
  parSCI[1] = -1;

  xpos  = 0.0;
  ypos  = 0.0;
  zpos  = 0.0;
  for (iplan = 0; iplan < kNplan; iplan++) {

    // The aluminum of the cross bars
    parSCB[0] = fCwidth[iplan]/2.0 + kSRLdst/2.0;
    sprintf(cTagV,"USF%01d",iplan);
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB);

    // The empty regions in the cross bars
    Float_t thkSCB = kSCBthk;
    if (iplan < 2) {
      thkSCB *= 1.5;
    }
    parSCI[2] = parSCB[2] - thkSCB;
    parSCI[0] = parSCB[0]/4.0 - kSCBthk;
    sprintf(cTagV,"USI%01d",iplan);
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parSCI,kNparSCI);

    sprintf(cTagV,"USI%01d",iplan);
    sprintf(cTagM,"USF%01d",iplan);
    ypos  = 0.0;
    zpos  = 0.0;
    xpos  =   parSCI[0] + thkSCB/2.0;
    gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
    xpos  = - parSCI[0] - thkSCB/2.0;
    gMC->Gspos(cTagV,2,cTagM,xpos,ypos,zpos,0,"ONLY");
    xpos  =   3.0 * parSCI[0] + 1.5 * thkSCB;
    gMC->Gspos(cTagV,3,cTagM,xpos,ypos,zpos,0,"ONLY");
    xpos  = - 3.0 * parSCI[0] - 1.5 * thkSCB;
    gMC->Gspos(cTagV,4,cTagM,xpos,ypos,zpos,0,"ONLY");

    sprintf(cTagV,"USF%01d",iplan);
    xpos  = 0.0;
    zpos  = fgkVrocsm + fgkSMpltT + parSCB[2] - fgkSheight/2.0 
          + iplan * (fgkCH + fgkVspace);

    ypos  =   fgkSlength/2.0 - kSCBwid/2.0;
    gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY");

    ypos  =   fClength[iplan][2]/2.0 + fClength[iplan][1];
    gMC->Gspos(cTagV,2,"UTI1", xpos,ypos,zpos,0,"ONLY");

    ypos  =   fClength[iplan][2]/2.0;
    gMC->Gspos(cTagV,3,"UTI1", xpos,ypos,zpos,0,"ONLY");

    ypos  = - fClength[iplan][2]/2.0;
    gMC->Gspos(cTagV,4,"UTI1", xpos,ypos,zpos,0,"ONLY");

    ypos  = - fClength[iplan][2]/2.0 - fClength[iplan][1];
    gMC->Gspos(cTagV,5,"UTI1", xpos,ypos,zpos,0,"ONLY");

    ypos  = - fgkSlength/2.0 + kSCBwid/2.0;
    gMC->Gspos(cTagV,6,"UTI1", xpos,ypos,zpos,0,"ONLY");

  }

  //
  // The horizontal connections between the cross bars
  //

  const Int_t   kNparSCH = 3;
  Float_t parSCH[kNparSCH];

  for (iplan = 1; iplan < kNplan-1; iplan++) {

    parSCH[0] = fCwidth[iplan]/2.0;
    parSCH[1] = (fClength[iplan+1][2]/2.0 + fClength[iplan+1][1]
               - fClength[iplan  ][2]/2.0 - fClength[iplan  ][1])/2.0;
    parSCH[2] = kSCHhgt/2.0;

    sprintf(cTagV,"USH%01d",iplan);
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCH,kNparSCH);
    xpos  = 0.0;
    ypos  = fClength[iplan][2]/2.0 + fClength[iplan][1] + parSCH[1];
    zpos  = fgkVrocsm + fgkSMpltT - kSCHhgt/2.0 - fgkSheight/2.0 
          + (iplan+1) * (fgkCH + fgkVspace);
    gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY");
    ypos  = -ypos;
    gMC->Gspos(cTagV,2,"UTI1", xpos,ypos,zpos,0,"ONLY");

  }

  //
  // The long corner ledges
  //

  const Int_t   kNparSCL  =  3;
  Float_t parSCL[kNparSCL];
  const Int_t   kNparSCLb = 11;
  Float_t parSCLb[kNparSCLb];

  // Upper ledges 
  // Thickness of the corner ledges
  const Float_t kSCLthkUa  =  0.6; 
  const Float_t kSCLthkUb  =  0.6; 
  // Width of the corner ledges
  const Float_t kSCLwidUa  =  3.2;
  const Float_t kSCLwidUb  =  4.8;
  // Position of the corner ledges
  const Float_t kSCLposxUa = 0.7;
  const Float_t kSCLposxUb = 3.3;
  const Float_t kSCLposzUa = 1.6;
  const Float_t kSCLposzUb = 0.3;
  // Vertical
  parSCL[0]  = kSCLthkUa /2.0;
  parSCL[1]  = fgkSlength/2.0;
  parSCL[2]  = kSCLwidUa /2.0;
  gMC->Gsvolu("USL1","BOX ",idtmed[1301-1],parSCL,kNparSCL);
  xpos  =   fgkSwidth2/2.0 - fgkSMpltT - kSCLposxUa;
  ypos  =   0.0;
  zpos  =   fgkSheight/2.0 - fgkSMpltT - kSCLposzUa;
  gMC->Gspos("USL1",1,"UTI1", xpos,ypos,zpos,matrix[0],"ONLY");
  xpos  = -xpos;
  gMC->Gspos("USL1",2,"UTI1", xpos,ypos,zpos,matrix[1],"ONLY");
  // Horizontal
  parSCL[0]  = kSCLwidUb /2.0;
  parSCL[1]  = fgkSlength/2.0;
  parSCL[2]  = kSCLthkUb /2.0;
  gMC->Gsvolu("USL2","BOX ",idtmed[1301-1],parSCL,kNparSCL);
  xpos  =   fgkSwidth2/2.0 - fgkSMpltT - kSCLposxUb;
  ypos  =   0.0;
  zpos  =   fgkSheight/2.0 - fgkSMpltT - kSCLposzUb; 
  gMC->Gspos("USL2",1,"UTI1", xpos,ypos,zpos,        0,"ONLY");
  xpos  = -xpos;
  gMC->Gspos("USL2",2,"UTI1", xpos,ypos,zpos,        0,"ONLY");

  // Lower ledges 
  // Thickness of the corner ledges
  const Float_t kSCLthkLa  =  2.464; 
  const Float_t kSCLthkLb  =  1.0; 
  // Width of the corner ledges
  const Float_t kSCLwidLa  =  8.5;
  const Float_t kSCLwidLb  =  3.3;
  // Position of the corner ledges
  const Float_t kSCLposxLa =  0.0;
  const Float_t kSCLposxLb =  2.6;
  const Float_t kSCLposzLa = -4.25;
  const Float_t kSCLposzLb = -0.5;
  // Vertical
  // Trapezoidal shape
  parSCLb[ 0] = fgkSlength/2.0;
  parSCLb[ 1] = 0.0;
  parSCLb[ 2] = 0.0;
  parSCLb[ 3] = kSCLwidLa /2.0;
  parSCLb[ 4] = kSCLthkLb /2.0;
  parSCLb[ 5] = kSCLthkLa /2.0;
  parSCLb[ 6] = 5.0;
  parSCLb[ 7] = kSCLwidLa /2.0;
  parSCLb[ 8] = kSCLthkLb /2.0;
  parSCLb[ 9] = kSCLthkLa /2.0;
  parSCLb[10] = 5.0;
  gMC->Gsvolu("USL3","TRAP",idtmed[1301-1],parSCLb,kNparSCLb);
  xpos  =   fgkSwidth1/2.0 - fgkSMpltT - kSCLposxLa;
  ypos  =   0.0;
  zpos  = - fgkSheight/2.0 + fgkSMpltT - kSCLposzLa;
  gMC->Gspos("USL3",1,"UTI1", xpos,ypos,zpos,matrix[2],"ONLY");
  xpos  = -xpos;
  gMC->Gspos("USL3",2,"UTI1", xpos,ypos,zpos,matrix[3],"ONLY");
  // Horizontal
  parSCL[0]  = kSCLwidLb /2.0;
  parSCL[1]  = fgkSlength/2.0;
  parSCL[2]  = kSCLthkLb /2.0;
  gMC->Gsvolu("USL4","BOX ",idtmed[1301-1],parSCL,kNparSCL);
  xpos  =   fgkSwidth1/2.0 - fgkSMpltT - kSCLposxLb;
  ypos  =   0.0;
  zpos  = - fgkSheight/2.0 + fgkSMpltT - kSCLposzLb;
  gMC->Gspos("USL4",1,"UTI1", xpos,ypos,zpos,        0,"ONLY");
  xpos  = -xpos;
  gMC->Gspos("USL4",2,"UTI1", xpos,ypos,zpos,        0,"ONLY");

}

//_____________________________________________________________________________
void AliTRDgeometry::CreateServices(Int_t *idtmed)
{
  //
  // Create the geometry of the services
  //
  // Names of the TRD services volumina
  //
  //        UTCL    Cooling arterias (Al)
  //        UTCW    Cooling arterias (Water)
  //        UUxx    Volumes for the services at the chambers (Air)
  //        UTPW    Power bars       (Cu)
  //        UTCP    Cooling pipes    (Fe)
  //        UTCH    Cooling pipes    (Water)
  //        UTPL    Power lines      (Cu)
  //        UMCM    Readout MCMs     (G10/Cu/Si)
  //

  Int_t   iplan = 0;
  Int_t   icham = 0;

  Float_t xpos  = 0.0;
  Float_t ypos  = 0.0;
  Float_t zpos  = 0.0;

  Char_t  cTagV[5];

  // The rotation matrices
  const Int_t kNmatrix = 4;
  Int_t   matrix[kNmatrix];
  gMC->Matrix(matrix[0], 100.0,   0.0,  90.0,  90.0,  10.0,   0.0);
  gMC->Matrix(matrix[1],  80.0,   0.0,  90.0,  90.0,  10.0, 180.0);
  gMC->Matrix(matrix[2],   0.0,   0.0,  90.0,  90.0,  90.0,   0.0);
  gMC->Matrix(matrix[3], 180.0,   0.0,  90.0,  90.0,  90.0, 180.0);

  AliTRDCommonParam *commonParam = AliTRDCommonParam::Instance();
  if (!commonParam) {
    AliError("Could not get common parameters\n");
    return;
  }
    
  //
  // The cooling arterias
  //

  // Width of the cooling arterias
  const Float_t kCOLwid  =  0.8; 
  // Height of the cooling arterias
  const Float_t kCOLhgt  =  6.5;
  // Positioning of the cooling 
  const Float_t kCOLposx =  1.8;
  const Float_t kCOLposz = -0.1;
  // Thickness of the walls of the cooling arterias
  const Float_t kCOLthk  =  0.1;
  const Int_t   kNparCOL =  3;
  Float_t parCOL[kNparCOL];
  parCOL[0]  = kCOLwid   /2.0;
  parCOL[1]  = fgkSlength/2.0;
  parCOL[2]  = kCOLhgt   /2.0;
  gMC->Gsvolu("UTCL","BOX ",idtmed[1308-1],parCOL,kNparCOL);
  parCOL[0] -= kCOLthk;
  parCOL[1]  = fgkSlength/2.0;
  parCOL[2] -= kCOLthk;
  gMC->Gsvolu("UTCW","BOX ",idtmed[1314-1],parCOL,kNparCOL);

  xpos  = 0.0;
  ypos  = 0.0;
  zpos  = 0.0;
  gMC->Gspos("UTCW",1,"UTCL", xpos,ypos,zpos,0,"ONLY");

  for (iplan = 1; iplan < kNplan; iplan++) { 

    xpos  = fCwidth[iplan]/2.0 + kCOLwid/2.0 + kCOLposx;
    ypos  = 0.0;
    zpos  = fgkVrocsm + fgkSMpltT + kCOLhgt/2.0 - fgkSheight/2.0 + kCOLposz 
          + iplan * (fgkCH + fgkVspace);
    gMC->Gspos("UTCL",iplan       ,"UTI1", xpos,ypos,zpos,matrix[0],"ONLY");
    gMC->Gspos("UTCL",iplan+kNplan,"UTI1",-xpos,ypos,zpos,matrix[1],"ONLY");

  }

  // The upper most layer (reaching into TOF acceptance)
  xpos  = fCwidth[5]/2.0 - kCOLhgt/2.0 - 1.3;
  ypos  = 0.0;
  zpos  = fgkSheight/2.0 - fgkSMpltT - 0.4 - kCOLwid/2.0; 
  gMC->Gspos("UTCL",6       ,"UTI1", xpos,ypos,zpos,matrix[3],"ONLY");
  gMC->Gspos("UTCL",6+kNplan,"UTI1",-xpos,ypos,zpos,matrix[3],"ONLY");

  //
  // The power bars
  //

  const Float_t kPWRwid  =  0.6;
  const Float_t kPWRhgt  =  5.0;
  const Float_t kPWRposx =  1.4;
  const Float_t kPWRposz =  1.9;
  const Int_t   kNparPWR =  3;
  Float_t parPWR[kNparPWR];
  parPWR[0] = kPWRwid   /2.0;
  parPWR[1] = fgkSlength/2.0;
  parPWR[2] = kPWRhgt   /2.0;
  gMC->Gsvolu("UTPW","BOX ",idtmed[1325-1],parPWR,kNparPWR);
  
  for (iplan = 1; iplan < kNplan; iplan++) { 
    
    xpos  = fCwidth[iplan]/2.0 + kPWRwid/2.0 + kPWRposx;
    ypos  = 0.0;
    zpos  = fgkVrocsm + fgkSMpltT + kPWRhgt/2.0 - fgkSheight/2.0 + kPWRposz 
          + iplan * (fgkCH + fgkVspace);
    gMC->Gspos("UTPW",iplan       ,"UTI1", xpos,ypos,zpos,matrix[0],"ONLY");
    gMC->Gspos("UTPW",iplan+kNplan,"UTI1",-xpos,ypos,zpos,matrix[1],"ONLY");

  }

  // The upper most layer (reaching into TOF acceptance)
  xpos  = fCwidth[5]/2.0 + kPWRhgt/2.0 - 1.3;
  ypos  = 0.0;
  zpos  = fgkSheight/2.0 - fgkSMpltT - 0.6 - kPWRwid/2.0; 
  gMC->Gspos("UTPW",6       ,"UTI1", xpos,ypos,zpos,matrix[3],"ONLY");
  gMC->Gspos("UTPW",6+kNplan,"UTI1",-xpos,ypos,zpos,matrix[3],"ONLY");

  //
  // The volumes for the services at the chambers
  //

  const Int_t kNparServ = 3;
  Float_t parServ[kNparServ];

  for (icham = 0; icham < kNcham; icham++) {
    for (iplan = 0; iplan < kNplan; iplan++) {

      Int_t iDet = GetDetectorSec(iplan,icham);

      sprintf(cTagV,"UU%02d",iDet);
      parServ[0] = fCwidth[iplan]        /2.0;
      parServ[1] = fClength[iplan][icham]/2.0 - fgkHspace/2.0;
      parServ[2] = fgkVspace             /2.0 - 0.742/2.0; 
      fChamberUUboxd[iDet][0] = parServ[0];
      fChamberUUboxd[iDet][1] = parServ[1];
      fChamberUUboxd[iDet][2] = parServ[2];
      gMC->Gsvolu(cTagV,"BOX",idtmed[1302-1],parServ,kNparServ);

      xpos  = 0.0;
      ypos  = - fClength[iplan][0] - fClength[iplan][1] - fClength[iplan][2]/2.0;
      for (Int_t ic = 0; ic < icham; ic++) {
        ypos += fClength[iplan][ic];        
      }
      ypos += fClength[iplan][icham]/2.0;
      zpos  = fgkVrocsm + fgkSMpltT + fgkCH + fgkVspace/2.0 - fgkSheight/2.0 
            + iplan * (fgkCH + fgkVspace);
      zpos -= 0.742/2.0;
      fChamberUUorig[iDet][0] = xpos;
      fChamberUUorig[iDet][1] = ypos;
      fChamberUUorig[iDet][2] = zpos;

    }
  }

  //
  // The cooling pipes inside the service volumes
  //

  const Int_t kNparTube = 3;
  Float_t parTube[kNparTube];
  // The cooling pipes
  parTube[0] = 0.0;
  parTube[1] = 0.0;
  parTube[2] = 0.0;
  gMC->Gsvolu("UTCP","TUBE",idtmed[1324-1],parTube,0);
  // The cooling water
  parTube[0] =  0.0;
  parTube[1] =  0.2/2.0;
  parTube[2] = -1.;
  gMC->Gsvolu("UTCH","TUBE",idtmed[1314-1],parTube,kNparTube);
  // Water inside the cooling pipe
  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  gMC->Gspos("UTCH",1,"UTCP",xpos,ypos,zpos,0,"ONLY");

  // Position the cooling pipes in the mother volume
  const Int_t kNpar = 3;
  Float_t par[kNpar];
  for (icham = 0; icham < kNcham;   icham++) {
    for (iplan = 0; iplan < kNplan; iplan++) {
      Int_t   iDet    = GetDetectorSec(iplan,icham);
      Int_t   iCopy   = GetDetector(iplan,icham,0) * 100;
      Int_t   nMCMrow = commonParam->GetRowMax(iplan,icham,0);
      Float_t ySize   = (GetChamberLength(iplan,icham) - 2.0*fgkRpadW) 
                      / ((Float_t) nMCMrow);
      sprintf(cTagV,"UU%02d",iDet);
      for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
        xpos   = 0.0;
        ypos   = (0.5 + iMCMrow) * ySize - 1.9 
               - fClength[iplan][icham]/2.0 + fgkHspace/2.0;
        zpos   = 0.0 + 0.742/2.0;                 
        par[0] = 0.0;
        par[1] = 0.3/2.0; // Thickness of the cooling pipes
        par[2] = fCwidth[iplan]/2.0;
        gMC->Gsposp("UTCP",iCopy+iMCMrow,cTagV,xpos,ypos,zpos
                          ,matrix[2],"ONLY",par,kNpar);
      }
    }
  }

  //
  // The power lines
  //

  // The copper power lines
  parTube[0] = 0.0;
  parTube[1] = 0.0;
  parTube[2] = 0.0;
  gMC->Gsvolu("UTPL","TUBE",idtmed[1305-1],parTube,0);

  // Position the power lines in the mother volume
  for (icham = 0; icham < kNcham;   icham++) {
    for (iplan = 0; iplan < kNplan; iplan++) {
      Int_t   iDet    = GetDetectorSec(iplan,icham);
      Int_t   iCopy   = GetDetector(iplan,icham,0) * 100;
      Int_t   nMCMrow = commonParam->GetRowMax(iplan,icham,0);
      Float_t ySize   = (GetChamberLength(iplan,icham) - 2.0*fgkRpadW) 
                      / ((Float_t) nMCMrow);
      sprintf(cTagV,"UU%02d",iDet);
      for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
        xpos   = 0.0;
        ypos   = (0.5 + iMCMrow) * ySize - 1.0 
               - fClength[iplan][icham]/2.0 + fgkHspace/2.0;
        zpos   = -0.4 + 0.742/2.0;
        par[0] = 0.0;
        par[1] = 0.2/2.0; // Thickness of the power lines
        par[2] = fCwidth[iplan]/2.0;
        gMC->Gsposp("UTPL",iCopy+iMCMrow,cTagV,xpos,ypos,zpos
                          ,matrix[2],"ONLY",par,kNpar);
      }
    }
  }

  //
  // The MCMs
  //

  const Float_t kMCMx    = 3.0;
  const Float_t kMCMy    = 3.0;
  const Float_t kMCMz    = 0.3;

  const Float_t kMCMpcTh = 0.1;
  const Float_t kMCMcuTh = 0.0215;
  const Float_t kMCMsiTh = 0.003;
  const Float_t kMCMcoTh = 0.1549;

  // The mother volume for the MCMs (air)
  const Int_t kNparMCM = 3;
  Float_t parMCM[kNparMCM];
  parMCM[0] = kMCMx   /2.0;
  parMCM[1] = kMCMy   /2.0;
  parMCM[2] = kMCMz   /2.0;
  gMC->Gsvolu("UMCM","BOX",idtmed[1302-1],parMCM,kNparMCM);

  // The MCM carrier G10 layer
  parMCM[0] = kMCMx   /2.0;
  parMCM[1] = kMCMy   /2.0;
  parMCM[2] = kMCMpcTh/2.0;
  gMC->Gsvolu("UMC1","BOX",idtmed[1319-1],parMCM,kNparMCM);
  // The MCM carrier Cu layer
  parMCM[0] = kMCMx   /2.0;
  parMCM[1] = kMCMy   /2.0;
  parMCM[2] = kMCMcuTh/2.0;
  gMC->Gsvolu("UMC2","BOX",idtmed[1318-1],parMCM,kNparMCM);
  // The silicon of the chips
  parMCM[0] = kMCMx   /2.0;
  parMCM[1] = kMCMy   /2.0;
  parMCM[2] = kMCMsiTh/2.0;
  gMC->Gsvolu("UMC3","BOX",idtmed[1320-1],parMCM,kNparMCM);
  // The aluminum of the cooling plates
  parMCM[0] = kMCMx   /2.0;
  parMCM[1] = kMCMy   /2.0;
  parMCM[2] = kMCMcoTh/2.0;
  gMC->Gsvolu("UMC4","BOX",idtmed[1324-1],parMCM,kNparMCM);

  // Put the MCM material inside the MCM mother volume
  xpos  =  0.0;
  ypos  =  0.0;
  zpos  = -kMCMz   /2.0 + kMCMpcTh/2.0;
  gMC->Gspos("UMC1",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  kMCMpcTh/2.0 + kMCMcuTh/2.0;
  gMC->Gspos("UMC2",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  kMCMcuTh/2.0 + kMCMsiTh/2.0;
  gMC->Gspos("UMC3",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  kMCMsiTh/2.0 + kMCMcoTh/2.0;
  gMC->Gspos("UMC4",1,"UMCM",xpos,ypos,zpos,0,"ONLY");

  // Position the MCMs in the mother volume
  for (icham = 0; icham < kNcham;   icham++) {
    for (iplan = 0; iplan < kNplan; iplan++) {
      Int_t   iDet    = GetDetectorSec(iplan,icham);
      Int_t   iCopy   = GetDetector(iplan,icham,0) * 1000;
      Int_t   nMCMrow = commonParam->GetRowMax(iplan,icham,0);
      Float_t ySize   = (GetChamberLength(iplan,icham) - 2.0*fgkRpadW) 
                      / ((Float_t) nMCMrow);
      Int_t   nMCMcol = 8;
      Float_t xSize   = (GetChamberWidth(iplan)        - 2.0*fgkCpadW)
	              / ((Float_t) nMCMcol);
      sprintf(cTagV,"UU%02d",iDet);
      for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
        for (Int_t iMCMcol = 0; iMCMcol < nMCMcol; iMCMcol++) {
          xpos   = (0.5 + iMCMcol) * xSize + 1.0 
                 - fCwidth[iplan]/2.0;
          ypos   = (0.5 + iMCMrow) * ySize + 1.0 
                 - fClength[iplan][icham]/2.0 + fgkHspace/2.0;
          zpos   = -0.4 + 0.742/2.0;
          par[0] = 0.0;
          par[1] = 0.2/2.0; // Thickness of the power lines
          par[2] = fCwidth[iplan]/2.0;
          gMC->Gspos("UMCM",iCopy+iMCMrow*10+iMCMcol,cTagV
                           ,xpos,ypos,zpos,0,"ONLY");
	}
      }

    }
  }

}

//_____________________________________________________________________________
void AliTRDgeometry::GroupChamber(Int_t iplan, Int_t icham, Int_t *idtmed)
{
  //
  // Group volumes UA, UD, UF, UU in a single chamber (Air)
  // UA, UD, UF, UU are boxes
  // UT will be a box
  //

  const Int_t kNparCha = 3;

  Int_t iDet = GetDetectorSec(iplan,icham);

  Float_t xyzMin[3];
  Float_t xyzMax[3];
  Float_t xyzOrig[3];
  Float_t xyzBoxd[3];

  Char_t  cTagV[5];
  Char_t  cTagM[5];

  for (Int_t i = 0; i < 3; i++) {
    xyzMin[i] = +9999.0; 
    xyzMax[i] = -9999.0;
  }

  for (Int_t i = 0; i < 3; i++) {

    xyzMin[i] = TMath::Min(xyzMin[i],fChamberUAorig[iDet][i]-fChamberUAboxd[iDet][i]);
    xyzMax[i] = TMath::Max(xyzMax[i],fChamberUAorig[iDet][i]+fChamberUAboxd[iDet][i]);

    xyzMin[i] = TMath::Min(xyzMin[i],fChamberUDorig[iDet][i]-fChamberUDboxd[iDet][i]);
    xyzMax[i] = TMath::Max(xyzMax[i],fChamberUDorig[iDet][i]+fChamberUDboxd[iDet][i]);

    xyzMin[i] = TMath::Min(xyzMin[i],fChamberUForig[iDet][i]-fChamberUFboxd[iDet][i]);
    xyzMax[i] = TMath::Max(xyzMax[i],fChamberUForig[iDet][i]+fChamberUFboxd[iDet][i]);

    xyzMin[i] = TMath::Min(xyzMin[i],fChamberUUorig[iDet][i]-fChamberUUboxd[iDet][i]);
    xyzMax[i] = TMath::Max(xyzMax[i],fChamberUUorig[iDet][i]+fChamberUUboxd[iDet][i]);

    xyzOrig[i] = 0.5*(xyzMax[i]+xyzMin[i]);
    xyzBoxd[i] = 0.5*(xyzMax[i]-xyzMin[i]);

  }
  
  sprintf(cTagM,"UT%02d",iDet);
  gMC->Gsvolu(cTagM,"BOX ",idtmed[1302-1],xyzBoxd,kNparCha);

  sprintf(cTagV,"UA%02d",iDet);
  gMC->Gspos(cTagV,1,cTagM
	    ,fChamberUAorig[iDet][0]-xyzOrig[0]
	    ,fChamberUAorig[iDet][1]-xyzOrig[1]
	    ,fChamberUAorig[iDet][2]-xyzOrig[2]
	    ,0,"ONLY");

  sprintf(cTagV,"UZ%02d",iDet);
  gMC->Gspos(cTagV,1,cTagM
	    ,fChamberUAorig[iDet][0]-xyzOrig[0] + fChamberUAboxd[iDet][0] - fgkCroW/2.0
	    ,fChamberUAorig[iDet][1]-xyzOrig[1]
	    ,fChamberUAorig[iDet][2]-xyzOrig[2] + fgkCraH/2.0 + fgkCdrH/2.0 - fgkCalW/2.0
	    ,0,"ONLY");
  gMC->Gspos(cTagV,2,cTagM
	    ,fChamberUAorig[iDet][0]-xyzOrig[0] - fChamberUAboxd[iDet][0] + fgkCroW/2.0
	    ,fChamberUAorig[iDet][1]-xyzOrig[1]
	    ,fChamberUAorig[iDet][2]-xyzOrig[2] + fgkCraH/2.0 + fgkCdrH/2.0 - fgkCalW/2.0
	    ,0,"ONLY");

  sprintf(cTagV,"UD%02d",iDet);
  gMC->Gspos(cTagV,1,cTagM
	    ,fChamberUDorig[iDet][0]-xyzOrig[0]
	    ,fChamberUDorig[iDet][1]-xyzOrig[1]
	    ,fChamberUDorig[iDet][2]-xyzOrig[2]
	    ,0,"ONLY");

  sprintf(cTagV,"UF%02d",iDet);
  gMC->Gspos(cTagV,1,cTagM
	    ,fChamberUForig[iDet][0]-xyzOrig[0]
	    ,fChamberUForig[iDet][1]-xyzOrig[1]
	    ,fChamberUForig[iDet][2]-xyzOrig[2]
	    ,0,"ONLY");
  
  sprintf(cTagV,"UU%02d",iDet);
  gMC->Gspos(cTagV,1,cTagM
            ,fChamberUUorig[iDet][0]-xyzOrig[0]
            ,fChamberUUorig[iDet][1]-xyzOrig[1]
            ,fChamberUUorig[iDet][2]-xyzOrig[2]
            ,0,"ONLY");

  sprintf(cTagV,"UT%02d",iDet);
  gMC->Gspos(cTagV,1,"UTI1"
            ,xyzOrig[0]
            ,xyzOrig[1]
            ,xyzOrig[2]
            ,0,"ONLY");

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

//CL
//_____________________________________________________________________________
Int_t AliTRDgeometry::GetPadRow(Int_t irob, Int_t imcm) const
{

  // return on which row this mcm sits 

  return fgkMCMrow*(irob/2) + imcm/fgkMCMrow;

;
}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetPadCol(Int_t irob, Int_t imcm, Int_t iadc) const
{
  //
  // return which pad is connected to this adc channel. return -1 if it
  // is one of the not directly connected adc channels (0, 1 20)
  //

  if (iadc < 2 || iadc > 19 ) return -1;

  return (iadc-2) + (imcm%fgkMCMrow)*fgkPadmax + GetRobSide(irob)*fgkColmax/2;

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetMCM(Int_t irow, Int_t icol) const
{

  // return on which mcm this pad is

  if ( irow < 0 || icol < 0 || irow > fgkRowmaxC1 || icol > fgkColmax ) return -1;

  return (icol%(fgkColmax/2))/fgkPadmax + fgkMCMrow*(irow%fgkMCMrow);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetROB(Int_t irow, Int_t icol) const
{

  // return on which rob this pad is

  return (irow/fgkMCMrow)*2 + GetColSide(icol);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetRobSide(Int_t irob) const
{

  // return on which side this rob sits (A side = 0, B side = 1)

  if ( irob < 0 || irob >= fgkROBmaxC1 ) return -1;

  return irob%2;

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetColSide(Int_t icol) const
{

  // return on which side this column sits (A side = 0, B side = 1)

  if ( icol < 0 || icol >= fgkColmax ) return -1;

  return icol/(fgkColmax/2);

}

//_____________________________________________________________________________
AliTRDgeometry *AliTRDgeometry::GetGeometry(AliRunLoader *runLoader)
{
  //
  // Load the geometry from the galice file
  //

  if (!runLoader) {
    runLoader = AliRunLoader::GetRunLoader();
  }
  if (!runLoader) {
    AliErrorGeneral("AliTRDgeometry::GetGeometry","No run loader");
    return NULL;
  }

  TDirectory *saveDir = gDirectory;
  runLoader->CdGAFile();

  // Try from the galice.root file
  AliTRDgeometry *geom = (AliTRDgeometry *) gDirectory->Get("TRDgeometry");

  if (!geom) {
    // If it is not in the file, try to get it from the run loader 
    AliTRD *trd = (AliTRD *) runLoader->GetAliRun()->GetDetector("TRD");
    geom = trd->GetGeometry();
  }
  if (!geom) {
    AliErrorGeneral("AliTRDgeometry::GetGeometry","Geometry not found");
    return NULL;
  }

  saveDir->cd();
  return geom;

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::ReadGeoMatrices()
{
  //
  // Read geo matrices from current gGeoManager for each TRD sector
  //

  if (!gGeoManager) {
    return kFALSE;
  }

  fMatrixArray           = new TObjArray(kNdet); 
  fMatrixCorrectionArray = new TObjArray(kNdet);
  fMatrixGeo             = new TObjArray(kNdet);
  AliAlignObjAngles o;

  for (Int_t iLayer = AliAlignObj::kTRD1; iLayer <= AliAlignObj::kTRD6; iLayer++) {
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(iLayer); iModule++) {

      UShort_t     volid   = AliAlignObj::LayerToVolUID(iLayer,iModule);
      const char  *symname = AliAlignObj::SymName(volid);
      TGeoPNEntry *pne     = gGeoManager->GetAlignableEntry(symname);
      const char  *path    = symname;
      if (pne) {
        path = pne->GetTitle();
      }
      if (!gGeoManager->cd(path)) {
        return kFALSE;
      }
      TGeoHMatrix *m         = gGeoManager->GetCurrentMatrix();
      Int_t        iLayerTRD = iLayer - AliAlignObj::kTRD1;
      Int_t        isector   = Nsect() - 1 - (iModule/Ncham());
      Int_t        ichamber  = Ncham() - 1 - (iModule%Ncham());
      Int_t        lid       = GetDetector(iLayerTRD,ichamber,isector);    

      //
      // Local geo system z-x-y  to x-y--z 
      //
      fMatrixGeo->AddAt(new TGeoHMatrix(*m),lid);
      
      TGeoRotation mchange; 
      mchange.RotateY(90); 
      mchange.RotateX(90);

      TGeoHMatrix gMatrix(mchange.Inverse());
      gMatrix.MultiplyLeft(m);
      fMatrixArray->AddAt(new TGeoHMatrix(gMatrix),lid); 

      //
      // Cluster transformation matrix
      //
      TGeoHMatrix  rotMatrix(mchange.Inverse());
      rotMatrix.MultiplyLeft(m);
      Double_t sectorAngle = 20.0 * (isector % 18) + 10.0;
      TGeoHMatrix  rotSector;
      rotSector.RotateZ(sectorAngle);
      rotMatrix.MultiplyLeft(&rotSector);      

      fMatrixCorrectionArray->AddAt(new TGeoHMatrix(rotMatrix),lid);       

    }    
  }

  return kTRUE;

}

