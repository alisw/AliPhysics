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
#include <TVirtualMC.h>
#include <TMath.h>

#include "AliLog.h"
#include "AliAlignObjParams.h"

#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"

ClassImp(AliTRDgeometry)

//_____________________________________________________________________________

  //
  // The geometry constants
  //
  const Int_t    AliTRDgeometry::fgkNsector   = kNsector;
  const Int_t    AliTRDgeometry::fgkNlayer    = kNlayer;
  const Int_t    AliTRDgeometry::fgkNstack    = kNstack;
  const Int_t    AliTRDgeometry::fgkNdet      = kNdet;

  //
  // Dimensions of the detector
  //

  // Total length of the TRD mother volume
  const Float_t  AliTRDgeometry::fgkTlength   = 751.0;

  // Parameter of the super module mother volumes 
  const Float_t  AliTRDgeometry::fgkSheight   =  77.9; 
  const Float_t  AliTRDgeometry::fgkSwidth1   =  94.881; 
  const Float_t  AliTRDgeometry::fgkSwidth2   = 122.353;
  const Float_t  AliTRDgeometry::fgkSlength   = 702.0;

  // Length of the additional space in front of the supermodule
  // used for services
  const Float_t  AliTRDgeometry::fgkFlength   = (AliTRDgeometry::fgkTlength
                                               - AliTRDgeometry::fgkSlength) / 2.0;

  // The super module side plates
  const Float_t  AliTRDgeometry::fgkSMpltT    =   0.2;

  // Vertical spacing of the chambers
  const Float_t  AliTRDgeometry::fgkVspace    =   1.784;
  // Horizontal spacing of the chambers
  const Float_t  AliTRDgeometry::fgkHspace    =   2.0;
  // Radial distance of the first ROC to the outer plates of the SM
  const Float_t  AliTRDgeometry::fgkVrocsm    =   1.2;

  // Height of different chamber parts
  // Radiator
  const Float_t  AliTRDgeometry::fgkCraH      =   4.8; 
  // Drift region
  const Float_t  AliTRDgeometry::fgkCdrH      =   3.0;
  // Amplification region
  const Float_t  AliTRDgeometry::fgkCamH      =   0.7;
  // Readout
  const Float_t  AliTRDgeometry::fgkCroH      =   2.316;
  // Additional width of the readout chamber frames
  const Float_t  AliTRDgeometry::fgkCroW      =   0.9;
  // Services on top of ROC
  const Float_t  AliTRDgeometry::fgkCsvH      = AliTRDgeometry::fgkVspace 
                                              -   0.742;
  // Total height (w/o services)
  const Float_t  AliTRDgeometry::fgkCH        = AliTRDgeometry::fgkCraH
                                              + AliTRDgeometry::fgkCdrH
                                              + AliTRDgeometry::fgkCamH
                                              + AliTRDgeometry::fgkCroH;  
  // Total height (with services)

  const Float_t  AliTRDgeometry::fgkCHsv      = AliTRDgeometry::fgkCH 
                                              + AliTRDgeometry::fgkCsvH;

  // Distance of anode wire plane relative to middle of alignable volume
  const Float_t  AliTRDgeometry::fgkAnodePos  = AliTRDgeometry::fgkCraH 
                                              + AliTRDgeometry::fgkCdrH 
                                              + AliTRDgeometry::fgkCamH/2.0
                                              - AliTRDgeometry::fgkCHsv/2.0;

  // Thicknesses of different parts of the chamber frame
  // Lower aluminum frame
  const Float_t  AliTRDgeometry::fgkCalT      =   0.4;
  // Lower Wacosit frame sides
  const Float_t  AliTRDgeometry::fgkCclsT     =   0.21;
  // Lower Wacosit frame front
  const Float_t  AliTRDgeometry::fgkCclfT     =   1.0;
  // Thickness of glue around radiator
  const Float_t  AliTRDgeometry::fgkCglT      =   0.25;
  // Upper Wacosit frame around amplification region
  const Float_t  AliTRDgeometry::fgkCcuTa     =   1.0;
  const Float_t  AliTRDgeometry::fgkCcuTb     =   0.8;
  // Al frame of back panel
  const Float_t  AliTRDgeometry::fgkCauT      =   1.5;
  // Additional Al ledge at the lower chamber frame
  // Actually the dimensions are not realistic, but 
  // modified in order to allow to mis-alignment. 
  // The amount of material is, however, correct 
  const Float_t  AliTRDgeometry::fgkCalW      =   2.5;
  const Float_t  AliTRDgeometry::fgkCalH      =   0.4;
  const Float_t  AliTRDgeometry::fgkCalWmod   =   0.4;
  const Float_t  AliTRDgeometry::fgkCalHmod   =   2.5;
  // Additional Wacosit ledge at the lower chamber frame
  const Float_t  AliTRDgeometry::fgkCwsW      =   1.2;
  const Float_t  AliTRDgeometry::fgkCwsH      =   0.3;

  // Difference of outer chamber width and pad plane width
  const Float_t  AliTRDgeometry::fgkCpadW     =   0.0;
  const Float_t  AliTRDgeometry::fgkRpadW     =   1.0;

  //
  // Thickness of the the material layers
  //
  const Float_t  AliTRDgeometry::fgkDrThick   = AliTRDgeometry::fgkCdrH;    
  const Float_t  AliTRDgeometry::fgkAmThick   = AliTRDgeometry::fgkCamH;
  const Float_t  AliTRDgeometry::fgkXeThick   = AliTRDgeometry::fgkDrThick
                                              + AliTRDgeometry::fgkAmThick;
  const Float_t  AliTRDgeometry::fgkWrThick   = 0.00011;

  const Float_t  AliTRDgeometry::fgkRMyThick  = 0.0015;
  const Float_t  AliTRDgeometry::fgkRCbThick  = 0.0055;
  const Float_t  AliTRDgeometry::fgkRGlThick  = 0.0065;
  const Float_t  AliTRDgeometry::fgkRRhThick  = 0.8;
  const Float_t  AliTRDgeometry::fgkRFbThick  = fgkCraH - 2.0 * (fgkRMyThick 
                                                               + fgkRCbThick 
                                                               + fgkRRhThick);

  const Float_t  AliTRDgeometry::fgkPPdThick  = 0.0025; 
  const Float_t  AliTRDgeometry::fgkPPpThick  = 0.0356; 
  const Float_t  AliTRDgeometry::fgkPGlThick  = 0.1428;
  const Float_t  AliTRDgeometry::fgkPCbThick  = 0.019;
  const Float_t  AliTRDgeometry::fgkPPcThick  = 0.0486;
  const Float_t  AliTRDgeometry::fgkPRbThick  = 0.0057;
  const Float_t  AliTRDgeometry::fgkPElThick  = 0.0029;
  const Float_t  AliTRDgeometry::fgkPHcThick  = fgkCroH - fgkPPdThick 
                                                        - fgkPPpThick
                                                        - fgkPGlThick 
                                                        - fgkPCbThick * 2.0
                                                        - fgkPPcThick
                                                        - fgkPRbThick
                                                        - fgkPElThick;

  //
  // Position of the material layers
  //
  const Float_t  AliTRDgeometry::fgkDrZpos    =  2.4;
  const Float_t  AliTRDgeometry::fgkAmZpos    =  0.0;
  const Float_t  AliTRDgeometry::fgkWrZposA   =  0.0;
  const Float_t  AliTRDgeometry::fgkWrZposB   = -fgkAmThick/2.0 + 0.001;
  const Float_t  AliTRDgeometry::fgkCalZpos   =  0.3;

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
  const Float_t  AliTRDgeometry::fgkTime0[6]  = { static_cast<Float_t>(fgkTime0Base + 0 * (Cheight() + Cspace())) 
                                                , static_cast<Float_t>(fgkTime0Base + 1 * (Cheight() + Cspace())) 
                                                , static_cast<Float_t>(fgkTime0Base + 2 * (Cheight() + Cspace())) 
                                                , static_cast<Float_t>(fgkTime0Base + 3 * (Cheight() + Cspace())) 
                                                , static_cast<Float_t>(fgkTime0Base + 4 * (Cheight() + Cspace())) 
                                                , static_cast<Float_t>(fgkTime0Base + 5 * (Cheight() + Cspace()))};

  const Double_t AliTRDgeometry::fgkXtrdBeg   = 288.43; // Values depend on position of TRD
  const Double_t AliTRDgeometry::fgkXtrdEnd   = 366.33; // mother volume inside space frame !!!

  // The outer width of the chambers
  const Float_t AliTRDgeometry::fgkCwidth[kNlayer] = {  90.4,  94.8,  99.3, 103.7, 108.1, 112.6 };
  
  // The outer lengths of the chambers
  // Includes the spacings between the chambers!
  const Float_t AliTRDgeometry::fgkClength[kNlayer][kNstack] = { { 124.0, 124.0, 110.0, 124.0, 124.0 }
							       , { 124.0, 124.0, 110.0, 124.0, 124.0 }
							       , { 131.0, 131.0, 110.0, 131.0, 131.0 }
							       , { 138.0, 138.0, 110.0, 138.0, 138.0 }
							       , { 145.0, 145.0, 110.0, 145.0, 145.0 }
							       , { 147.0, 147.0, 110.0, 147.0, 147.0 } };

        Char_t  AliTRDgeometry::fgSMstatus[kNsector]         = { 1, 1, 1, 1, 1, 1, 1, 1, 1
				       		               , 1, 1, 1, 1, 1, 1, 1, 1, 1 };

  TObjArray* AliTRDgeometry::fgClusterMatrixArray = NULL;

  TObjArray* AliTRDgeometry::fgPadPlaneArray = NULL;

//_____________________________________________________________________________
AliTRDgeometry::AliTRDgeometry()
{
  //
  // AliTRDgeometry default constructor
  //

}

//_____________________________________________________________________________
AliTRDgeometry::~AliTRDgeometry()
{
  //
  // AliTRDgeometry destructor
  //

}

//_____________________________________________________________________________
void AliTRDgeometry::CreatePadPlaneArray()
{
  //
  // Creates the array of AliTRDpadPlane objects
  //

  if (fgPadPlaneArray)
    return;

  static TObjArray padPlaneArray(fgkNlayer * fgkNstack);
  padPlaneArray.SetOwner(kTRUE);

  fgPadPlaneArray = &padPlaneArray;
  for (Int_t ilayer = 0; ilayer < fgkNlayer; ilayer++) {
    for (Int_t istack = 0; istack < fgkNstack; istack++) {
      Int_t ipp = GetDetectorSec(ilayer,istack);
      fgPadPlaneArray->AddAt(CreatePadPlane(ilayer,istack),ipp);
    }
  }

}

//_____________________________________________________________________________
AliTRDpadPlane *AliTRDgeometry::CreatePadPlane(Int_t ilayer, Int_t istack)
{
  //
  // Creates an AliTRDpadPlane object
  //

  AliTRDpadPlane *padPlane = new AliTRDpadPlane();

  padPlane->SetLayer(ilayer);
  padPlane->SetStack(istack);

  padPlane->SetRowSpacing(0.0);
  padPlane->SetColSpacing(0.0);

  padPlane->SetLengthRim(1.0);
  padPlane->SetWidthRim(0.5);

  padPlane->SetNcols(144);

  padPlane->SetAnodeWireOffset(0.25);

  //
  // The pad plane parameter
  //
  const Float_t kTiltAngle = 2.0;
  switch (ilayer) {
  case 0:
    if (istack == 2) {
      // L0C0 type
      padPlane->SetNrows(12);
      padPlane->SetLength(108.0);
      padPlane->SetLengthOPad(8.0);
      padPlane->SetLengthIPad(9.0);
    }
    else {
      // L0C1 type
      padPlane->SetNrows(16);
      padPlane->SetLength(122.0);
      padPlane->SetLengthOPad(7.5);
      padPlane->SetLengthIPad(7.5);
    }
    padPlane->SetWidth(92.2);
    padPlane->SetWidthOPad(0.515);
    padPlane->SetWidthIPad(0.635);
    padPlane->SetTiltingAngle(-kTiltAngle);
    break;
  case 1:
    if (istack == 2) {
      // L1C0 type
      padPlane->SetNrows(12);
      padPlane->SetLength(108.0);
      padPlane->SetLengthOPad(8.0);
      padPlane->SetLengthIPad(9.0);
    }
    else {
      // L1C1 type
      padPlane->SetNrows(16);
      padPlane->SetLength(122.0);
      padPlane->SetLengthOPad(7.5);
      padPlane->SetLengthIPad(7.5);
    }
    padPlane->SetWidth(96.6);
    padPlane->SetWidthOPad(0.585);
    padPlane->SetWidthIPad(0.665);
    padPlane->SetTiltingAngle(kTiltAngle);
    break;
  case 2:
    if (istack == 2) {
      // L2C0 type
      padPlane->SetNrows(12);
      padPlane->SetLength(108.0);
      padPlane->SetLengthOPad(8.0);
      padPlane->SetLengthIPad(9.0);
    }
    else {
      // L2C1 type
      padPlane->SetNrows(16);
      padPlane->SetLength(129.0);
      padPlane->SetLengthOPad(7.5);
      padPlane->SetLengthIPad(8.0);
    }
    padPlane->SetWidth(101.1);
    padPlane->SetWidthOPad(0.705);
    padPlane->SetWidthIPad(0.695);
    padPlane->SetTiltingAngle(-kTiltAngle);
    break;
  case 3:
    if (istack == 2) {
      // L3C0 type
      padPlane->SetNrows(12);
      padPlane->SetLength(108.0);
      padPlane->SetLengthOPad(8.0);
      padPlane->SetLengthIPad(9.0);
    }
    else {
      // L3C1 type
      padPlane->SetNrows(16);
      padPlane->SetLength(136.0);
      padPlane->SetLengthOPad(7.5);
      padPlane->SetLengthIPad(8.5);
    }
    padPlane->SetWidth(105.5);
    padPlane->SetWidthOPad(0.775);
    padPlane->SetWidthIPad(0.725);
    padPlane->SetTiltingAngle(kTiltAngle);
    break;
  case 4:
    if (istack == 2) {
      // L4C0 type
      padPlane->SetNrows(12);
      padPlane->SetLength(108.0);
      padPlane->SetLengthOPad(8.0);
    }
    else {
      // L4C1 type
      padPlane->SetNrows(16);
      padPlane->SetLength(143.0);
      padPlane->SetLengthOPad(7.5);
    }
    padPlane->SetWidth(109.9);
    padPlane->SetWidthOPad(0.845);
    padPlane->SetLengthIPad(9.0);
    padPlane->SetWidthIPad(0.755);
    padPlane->SetTiltingAngle(-kTiltAngle);
    break;
  case 5:
    if (istack == 2) {
      // L5C0 type
      padPlane->SetNrows(12);
      padPlane->SetLength(108.0);
      padPlane->SetLengthOPad(8.0);
    }
    else {
      // L5C1 type
      padPlane->SetNrows(16);
      padPlane->SetLength(145.0);
      padPlane->SetLengthOPad(8.5);
    }
    padPlane->SetWidth(114.4);
    padPlane->SetWidthOPad(0.965);
    padPlane->SetLengthIPad(9.0);
    padPlane->SetWidthIPad(0.785);
    padPlane->SetTiltingAngle(kTiltAngle);
    break;
  };

  //
  // The positions of the borders of the pads
  //
  // Row direction
  //
  Double_t row = fgkClength[ilayer][istack] / 2.0
               - fgkRpadW
               - padPlane->GetLengthRim();
  for (Int_t ir = 0; ir < padPlane->GetNrows(); ir++) {
    padPlane->SetPadRow(ir,row);
    row -= padPlane->GetRowSpacing();
    if (ir == 0) {
      row -= padPlane->GetLengthOPad();
    }
    else {
      row -= padPlane->GetLengthIPad();
    }
  }
  //
  // Column direction
  //
  Double_t col = - fgkCwidth[ilayer] / 2.0
                 - fgkCroW
                 + padPlane->GetWidthRim();
  for (Int_t ic = 0; ic < padPlane->GetNcols(); ic++) {
    padPlane->SetPadCol(ic,col);
    col += padPlane->GetColSpacing();
    if (ic == 0) {
      col += padPlane->GetWidthOPad();
    }
    else {
      col += padPlane->GetWidthIPad();
    }
  }
  // Calculate the offset to translate from the local ROC system into
  // the local supermodule system, which is used for clusters
  Double_t rowTmp = fgkClength[ilayer][0]
    	          + fgkClength[ilayer][1]
                  + fgkClength[ilayer][2] / 2.0;
  for (Int_t jstack = 0; jstack < istack; jstack++) {
    rowTmp -= fgkClength[ilayer][jstack];
  }
  padPlane->SetPadRowSMOffset(rowTmp - fgkClength[ilayer][istack]/2.0);

  return padPlane;

}

//_____________________________________________________________________________
void AliTRDgeometry::CreateGeometry(Int_t *idtmed)
{
  //
  // Create the TRD geometry
  //
  //
  // Names of the TRD volumina (xx = detector number):
  //
  //   Volume (Air) wrapping the readout chamber components
  //     UTxx    includes: UAxx, UDxx, UFxx, UUxx
  //
  //   Lower part of the readout chambers (drift volume + radiator)
  //     UAxx    Aluminum frames                (Al)
  //
  //   Upper part of the readout chambers (readout plane + fee)
  //     UDxx    Wacosit frames of amp. region  (Wacosit)
  //     UFxx    Aluminum frame of back panel   (Al)
  //
  //   Services on chambers (cooling, cables, MCMs, DCS boards, ...)
  //     UUxx    Volume containing the services (Air) 
  //
  //   Material layers inside sensitive area:
  //     Name    Description                     Mat.      Thick.   Dens.    Radl.    X/X_0
  //                                                        
  //     URMYxx  Mylar layers (x2)               Mylar     0.0015   1.39     28.5464  0.005%
  //     URCBxx  Carbon layer (x2)               Carbon    0.0055   1.75     24.2824  0.023%
  //     URGLxx  Glue on the carbon layers (x2)  Araldite  0.0065   1.12     37.0664  0.018%
  //     URRHxx  Rohacell layer (x2)             Rohacell  0.8      0.075    536.005  0.149%
  //     URFBxx  Fiber mat layer                 PP        3.186    0.068    649.727  0.490%
  //     
  //     UJxx    Drift region                    Xe/CO2    3.0      0.00495  1792.37  0.167%
  //     UKxx    Amplification region            Xe/CO2    0.7      0.00495  1792.37  0.039%
  //     UWxx    Wire planes (x2)                Copper    0.00011  8.96     1.43503  0.008%
  //
  //     UPPDxx  Copper of pad plane             Copper    0.0025   8.96     1.43503  0.174%
  //     UPPPxx  PCB of pad plane                G10       0.0356   2.0      14.9013  0.239%
  //     UPGLxx  Glue on pad planes              Araldite  0.0923   1.12     37.0664  0.249%
  //             + add. glue (ca. 600g)          Araldite  0.0505   1.12     37.0663  0.107%
  //     UPCBxx  Carbon fiber mats (x2)          Carbon    0.019    1.75     24.2824  0.078%
  //     UPHCxx  Honeycomb structure             Aramide   2.0299   0.032    1198.84  0.169%
  //     UPPCxx  PCB of readout board            G10       0.0486   2.0      14.9013  0.326%
  //     UPRDxx  Copper of readout board         Copper    0.0057   8.96     1.43503  0.404%
  //     UPELxx  Electronics + cables            Copper    0.0029   8.96     1.43503  0.202%
  //

  const Int_t kNparTrd = 4;
  const Int_t kNparCha = 3;

  Float_t xpos;
  Float_t ypos;
  Float_t zpos;

  Float_t parTrd[kNparTrd];
  Float_t parCha[kNparCha];

  const Int_t kTag = 100;
  Char_t  cTagV[kTag];
  Char_t  cTagM[kTag];

  // There are three TRD volumes for the supermodules in order to accomodate
  // the different arrangements in front of PHOS
  // UTR1: Default supermodule
  // UTR2: Supermodule in front of PHOS with double carbon cover
  // UTR3: As UTR2, but w/o middle stack
  // UTR4: Sector 17 with missing chamber L4S4
  //
  // The mother volume for one sector (Air), full length in z-direction
  // Provides material for side plates of super module
  parTrd[0] = fgkSwidth1/2.0;
  parTrd[1] = fgkSwidth2/2.0;
  parTrd[2] = fgkSlength/2.0;
  parTrd[3] = fgkSheight/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTR1","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  TVirtualMC::GetMC()->Gsvolu("UTR2","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  TVirtualMC::GetMC()->Gsvolu("UTR3","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  TVirtualMC::GetMC()->Gsvolu("UTR4","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  // The outer aluminum plates of the super module (Al)
  parTrd[0] = fgkSwidth1/2.0;
  parTrd[1] = fgkSwidth2/2.0;
  parTrd[2] = fgkSlength/2.0;
  parTrd[3] = fgkSheight/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTS1","TRD1",idtmed[1301-1],parTrd,kNparTrd);
  TVirtualMC::GetMC()->Gsvolu("UTS2","TRD1",idtmed[1301-1],parTrd,kNparTrd);
  TVirtualMC::GetMC()->Gsvolu("UTS3","TRD1",idtmed[1301-1],parTrd,kNparTrd);
  TVirtualMC::GetMC()->Gsvolu("UTS4","TRD1",idtmed[1301-1],parTrd,kNparTrd);
  // The inner part of the TRD mother volume for one sector (Air), 
  // full length in z-direction
  parTrd[0] = fgkSwidth1/2.0 - fgkSMpltT;
  parTrd[1] = fgkSwidth2/2.0 - fgkSMpltT;
  parTrd[2] = fgkSlength/2.0;
  parTrd[3] = fgkSheight/2.0 - fgkSMpltT;
  TVirtualMC::GetMC()->Gsvolu("UTI1","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  TVirtualMC::GetMC()->Gsvolu("UTI2","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  TVirtualMC::GetMC()->Gsvolu("UTI3","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  TVirtualMC::GetMC()->Gsvolu("UTI4","TRD1",idtmed[1302-1],parTrd,kNparTrd);

  // The inner part of the TRD mother volume for services in front
  // of the supermodules  (Air), 
  parTrd[0] = fgkSwidth1/2.0;
  parTrd[1] = fgkSwidth2/2.0;
  parTrd[2] = fgkFlength/2.0;
  parTrd[3] = fgkSheight/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTF1","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  TVirtualMC::GetMC()->Gsvolu("UTF2","TRD1",idtmed[1302-1],parTrd,kNparTrd);

  for (Int_t istack = 0; istack < kNstack; istack++) {
    for (Int_t ilayer = 0; ilayer < kNlayer; ilayer++) {  

      Int_t iDet = GetDetectorSec(ilayer,istack);

      // The lower part of the readout chambers (drift volume + radiator) 
      // The aluminum frames 
      snprintf(cTagV,kTag,"UA%02d",iDet);
      parCha[0] = fgkCwidth[ilayer]/2.0;
      parCha[1] = fgkClength[ilayer][istack]/2.0 - fgkHspace/2.0;
      parCha[2] = fgkCraH/2.0 + fgkCdrH/2.0;
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
      // The additional aluminum on the frames
      // This part has not the correct shape but is just supposed to
      // represent the missing material. The correct form of the L-shaped
      // profile would not fit into the alignable volume. 
      snprintf(cTagV,kTag,"UZ%02d",iDet);
      parCha[0] = fgkCalWmod/2.0;
      parCha[1] = fgkClength[ilayer][istack]/2.0 - fgkHspace/2.0;
      parCha[2] = fgkCalHmod/2.0;
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
      // The additional Wacosit on the frames
      snprintf(cTagV,kTag,"UP%02d",iDet);
      parCha[0] = fgkCwsW/2.0;
      parCha[1] = fgkClength[ilayer][istack]/2.0 - fgkHspace/2.0;
      parCha[2] = fgkCwsH/2.0;
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
      // The Wacosit frames 
      snprintf(cTagV,kTag,"UB%02d",iDet);
      parCha[0] = fgkCwidth[ilayer]/2.0 - fgkCalT; 
      parCha[1] = -1.0;
      parCha[2] = -1.0;
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
      // The glue around the radiator
      snprintf(cTagV,kTag,"UX%02d",iDet);
      parCha[0] = fgkCwidth[ilayer]/2.0 - fgkCalT - fgkCclsT; 
      parCha[1] = fgkClength[ilayer][istack]/2.0 - fgkHspace/2.0 - fgkCclfT;
      parCha[2] = fgkCraH/2.0;
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1311-1],parCha,kNparCha);
      // The inner part of radiator (air)
      snprintf(cTagV,kTag,"UC%02d",iDet);
      parCha[0] = fgkCwidth[ilayer]/2.0 - fgkCalT - fgkCclsT - fgkCglT; 
      parCha[1] = fgkClength[ilayer][istack]/2.0 - fgkHspace/2.0 - fgkCclfT - fgkCglT;
      parCha[2] = -1.0;
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);

      // The upper part of the readout chambers (amplification volume)
      // The Wacosit frames
      snprintf(cTagV,kTag,"UD%02d",iDet);
      parCha[0] = fgkCwidth[ilayer]/2.0 + fgkCroW;
      parCha[1] = fgkClength[ilayer][istack]/2.0 - fgkHspace/2.0;
      parCha[2] = fgkCamH/2.0;
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
      // The inner part of the Wacosit frame (air)
      snprintf(cTagV,kTag,"UE%02d",iDet);
      parCha[0] = fgkCwidth[ilayer]/2.0 + fgkCroW - fgkCcuTb; 
      parCha[1] = fgkClength[ilayer][istack]/2.0 - fgkHspace/2.0 - fgkCcuTa;
      parCha[2] = -1.;
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);

      // The back panel, including pad plane and readout boards
      // The aluminum frames
      snprintf(cTagV,kTag,"UF%02d",iDet);
      parCha[0] = fgkCwidth[ilayer]/2.0 + fgkCroW;
      parCha[1] = fgkClength[ilayer][istack]/2.0 - fgkHspace/2.0;
      parCha[2] = fgkCroH/2.0;
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
      // The inner part of the aluminum frames
      snprintf(cTagV,kTag,"UG%02d",iDet);
      parCha[0] = fgkCwidth[ilayer]/2.0 + fgkCroW - fgkCauT; 
      parCha[1] = fgkClength[ilayer][istack]/2.0 - fgkHspace/2.0 - fgkCauT;
      parCha[2] = -1.0;
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);

      //
      // The material layers inside the chambers
      //

      // Mylar layer (radiator)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkRMyThick/2.0;
      snprintf(cTagV,kTag,"URMY%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1327-1],parCha,kNparCha);
      // Carbon layer (radiator)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkRCbThick/2.0;
      snprintf(cTagV,kTag,"URCB%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1326-1],parCha,kNparCha);
      // Araldite layer (radiator)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkRGlThick/2.0;
      snprintf(cTagV,kTag,"URGL%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1311-1],parCha,kNparCha);
      // Rohacell layer (radiator)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkRRhThick/2.0;
      snprintf(cTagV,kTag,"URRH%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1315-1],parCha,kNparCha);
      // Fiber layer (radiator)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkRFbThick/2.0;
      snprintf(cTagV,kTag,"URFB%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1328-1],parCha,kNparCha);

      // Xe/Isobutane layer (drift volume) 
      parCha[0] = fgkCwidth[ilayer]/2.0 - fgkCalT - fgkCclsT;
      parCha[1] = fgkClength[ilayer][istack]/2.0 - fgkHspace/2.0 - fgkCclfT;
      parCha[2] = fgkDrThick/2.0;
      snprintf(cTagV,kTag,"UJ%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1309-1],parCha,kNparCha);

      // Xe/Isobutane layer (amplification volume)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkAmThick/2.0;
      snprintf(cTagV,kTag,"UK%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1309-1],parCha,kNparCha);  
      // Cu layer (wire plane)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkWrThick/2.0;
      snprintf(cTagV,kTag,"UW%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1303-1],parCha,kNparCha);

      // Cu layer (pad plane)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkPPdThick/2.0;
      snprintf(cTagV,kTag,"UPPD%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1305-1],parCha,kNparCha);
      // G10 layer (pad plane)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkPPpThick/2.0;
      snprintf(cTagV,kTag,"UPPP%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1313-1],parCha,kNparCha);
      // Araldite layer (glue)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkPGlThick/2.0;
      snprintf(cTagV,kTag,"UPGL%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1311-1],parCha,kNparCha);
      // Carbon layer (carbon fiber mats)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkPCbThick/2.0;
      snprintf(cTagV,kTag,"UPCB%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1326-1],parCha,kNparCha);
      // Aramide layer (honeycomb)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkPHcThick/2.0;
      snprintf(cTagV,kTag,"UPHC%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1310-1],parCha,kNparCha);
      // G10 layer (PCB readout board)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkPPcThick/2;
      snprintf(cTagV,kTag,"UPPC%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1313-1],parCha,kNparCha);
      // Cu layer (traces in readout board)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkPRbThick/2.0;
      snprintf(cTagV,kTag,"UPRB%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1306-1],parCha,kNparCha);
      // Cu layer (other material on in readout board, incl. screws)
      parCha[0] = -1.0;
      parCha[1] = -1.0;
      parCha[2] = fgkPElThick/2.0;
      snprintf(cTagV,kTag,"UPEL%02d",iDet);
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1304-1],parCha,kNparCha);

      //
      // Position the layers in the chambers
      //
      xpos = 0.0;
      ypos = 0.0;

      // Lower part
      // Mylar layers (radiator)
      zpos =  fgkRMyThick/2.0 - fgkCraH/2.0;
      snprintf(cTagV,kTag,"URMY%02d",iDet);
      snprintf(cTagM,kTag,"UC%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      zpos = -fgkRMyThick/2.0 + fgkCraH/2.0;
      snprintf(cTagV,kTag,"URMY%02d",iDet);
      snprintf(cTagM,kTag,"UC%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,2,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Carbon layers (radiator)
      zpos =  fgkRCbThick/2.0 + fgkRMyThick - fgkCraH/2.0;
      snprintf(cTagV,kTag,"URCB%02d",iDet);
      snprintf(cTagM,kTag,"UC%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      zpos = -fgkRCbThick/2.0 - fgkRMyThick + fgkCraH/2.0;
      snprintf(cTagV,kTag,"URCB%02d",iDet);
      snprintf(cTagM,kTag,"UC%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,2,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Carbon layers (radiator)
      zpos =  fgkRGlThick/2.0 + fgkRCbThick + fgkRMyThick - fgkCraH/2.0;
      snprintf(cTagV,kTag,"URGL%02d",iDet);
      snprintf(cTagM,kTag,"UC%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      zpos = -fgkRGlThick/2.0 - fgkRCbThick - fgkRMyThick + fgkCraH/2.0;
      snprintf(cTagV,kTag,"URGL%02d",iDet);
      snprintf(cTagM,kTag,"UC%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,2,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Rohacell layers (radiator)
      zpos =  fgkRRhThick/2.0 + fgkRGlThick + fgkRCbThick + fgkRMyThick - fgkCraH/2.0;
      snprintf(cTagV,kTag,"URRH%02d",iDet);
      snprintf(cTagM,kTag,"UC%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      zpos = -fgkRRhThick/2.0 - fgkRGlThick - fgkRCbThick - fgkRMyThick + fgkCraH/2.0;
      snprintf(cTagV,kTag,"URRH%02d",iDet);
      snprintf(cTagM,kTag,"UC%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,2,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Fiber layers (radiator)
      zpos =  0.0;
      snprintf(cTagV,kTag,"URFB%02d",iDet);
      snprintf(cTagM,kTag,"UC%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");

      // Xe/Isobutane layer (drift volume) 
      zpos = fgkDrZpos;
      snprintf(cTagV,kTag,"UJ%02d",iDet);
      snprintf(cTagM,kTag,"UB%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");

      // Upper part
      // Xe/Isobutane layer (amplification volume)
      zpos = fgkAmZpos;
      snprintf(cTagV,kTag,"UK%02d",iDet);
      snprintf(cTagM,kTag,"UE%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Cu layer (wire planes inside amplification volume)
      zpos = fgkWrZposA; 
      snprintf(cTagV,kTag,"UW%02d",iDet);
      snprintf(cTagM,kTag,"UK%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      zpos = fgkWrZposB; 
      snprintf(cTagV,kTag,"UW%02d",iDet);
      snprintf(cTagM,kTag,"UK%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,2,cTagM,xpos,ypos,zpos,0,"ONLY");

      // Back panel + pad plane + readout part
      // Cu layer (pad plane)
      zpos =  fgkPPdThick/2.0 - fgkCroH/2.0;
      snprintf(cTagV,kTag,"UPPD%02d",iDet);
      snprintf(cTagM,kTag,"UG%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // G10  layer (pad plane)
      zpos =  fgkPPpThick/2.0 + fgkPPdThick - fgkCroH/2.0;
      snprintf(cTagV,kTag,"UPPP%02d",iDet);
      snprintf(cTagM,kTag,"UG%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Araldite layer (glue)
      zpos =  fgkPGlThick/2.0 + fgkPPpThick + fgkPPdThick - fgkCroH/2.0;
      snprintf(cTagV,kTag,"UPGL%02d",iDet);
      snprintf(cTagM,kTag,"UG%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Carbon layers (carbon fiber mats)
      zpos =  fgkPCbThick/2.0 + fgkPGlThick + fgkPPpThick + fgkPPdThick - fgkCroH/2.0;
      snprintf(cTagV,kTag,"UPCB%02d",iDet);
      snprintf(cTagM,kTag,"UG%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      zpos = -fgkPCbThick/2.0 - fgkPPcThick - fgkPRbThick - fgkPElThick + fgkCroH/2.0;
      snprintf(cTagV,kTag,"UPCB%02d",iDet);
      snprintf(cTagM,kTag,"UG%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,2,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Aramide layer (honeycomb)
      zpos =  fgkPHcThick/2.0 + fgkPCbThick + fgkPGlThick + fgkPPpThick + fgkPPdThick - fgkCroH/2.0;
      snprintf(cTagV,kTag,"UPHC%02d",iDet);
      snprintf(cTagM,kTag,"UG%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // G10 layer (PCB readout board)
      zpos = -fgkPPcThick/2.0 - fgkPRbThick - fgkPElThick + fgkCroH/2.0;
      snprintf(cTagV,kTag,"UPPC%02d",iDet);
      snprintf(cTagM,kTag,"UG%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Cu layer (traces in readout board)
      zpos = -fgkPRbThick/2.0 - fgkPElThick + fgkCroH/2.0;
      snprintf(cTagV,kTag,"UPRB%02d",iDet);
      snprintf(cTagM,kTag,"UG%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Cu layer (other materials on readout board, incl. screws)
      zpos = -fgkPElThick/2.0 + fgkCroH/2.0;
      snprintf(cTagV,kTag,"UPEL%02d",iDet);
      snprintf(cTagM,kTag,"UG%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");

      // Position the inner volumes of the chambers in the frames
      xpos = 0.0;
      ypos = 0.0;

      // The inner part of the radiator (air)
      zpos = 0.0;
      snprintf(cTagV,kTag,"UC%02d",iDet);
      snprintf(cTagM,kTag,"UX%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // The glue around the radiator
      zpos = fgkCraH/2.0 - fgkCdrH/2.0 - fgkCraH/2.0;
      snprintf(cTagV,kTag,"UX%02d",iDet);
      snprintf(cTagM,kTag,"UB%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // The lower Wacosit frame inside the aluminum frame
      zpos = 0.0;
      snprintf(cTagV,kTag,"UB%02d",iDet);
      snprintf(cTagM,kTag,"UA%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");

      // The inside of the upper Wacosit frame
      zpos = 0.0;
      snprintf(cTagV,kTag,"UE%02d",iDet);
      snprintf(cTagM,kTag,"UD%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");

      // The inside of the upper aluminum frame
      zpos = 0.0;
      snprintf(cTagV,kTag,"UG%02d",iDet);
      snprintf(cTagM,kTag,"UF%02d",iDet);
      TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");      

    }
  }

  // Create the volumes of the super module frame
  CreateFrame(idtmed);

  // Create the volumes of the services
  CreateServices(idtmed);
  
  for (Int_t istack = 0; istack < kNstack; istack++) {
    for (Int_t ilayer = 0; ilayer < kNlayer; ilayer++) {  
      AssembleChamber(ilayer,istack);
    }
  }
  
  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  TVirtualMC::GetMC()->Gspos("UTI1",1,"UTS1",xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTI2",1,"UTS2",xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTI3",1,"UTS3",xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTI4",1,"UTS4",xpos,ypos,zpos,0,"ONLY");

  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  TVirtualMC::GetMC()->Gspos("UTS1",1,"UTR1",xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTS2",1,"UTR2",xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTS3",1,"UTR3",xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTS4",1,"UTR4",xpos,ypos,zpos,0,"ONLY");

  // Put the TRD volumes into the space frame mother volumes
  // if enabled via status flag
  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  for (Int_t isector = 0; isector < kNsector; isector++) {
    if (GetSMstatus(isector)) {
      snprintf(cTagV,kTag,"BTRD%d",isector);
      switch (isector) {
      case 17:
        // Missing L4S4 chamber
        TVirtualMC::GetMC()->Gspos("UTR4",1,cTagV,xpos,ypos,zpos,0,"ONLY");
        break;
      case 13:
      case 14:
      case 15:
        // Double carbon, w/o middle stack
        TVirtualMC::GetMC()->Gspos("UTR3",1,cTagV,xpos,ypos,zpos,0,"ONLY");
        break;
      case 11:
      case 12:
	// Double carbon, all stacks
        TVirtualMC::GetMC()->Gspos("UTR2",1,cTagV,xpos,ypos,zpos,0,"ONLY");
        break;
      default:
	// Standard supermodule
        TVirtualMC::GetMC()->Gspos("UTR1",1,cTagV,xpos,ypos,zpos,0,"ONLY");
      };
    }
  }

  // Put the TRD volumes into the space frame mother volumes
  // if enabled via status flag
  xpos = 0.0;
  ypos = 0.5*fgkSlength + 0.5*fgkFlength;
  zpos = 0.0;
  for (Int_t isector = 0; isector < kNsector; isector++) {
    if (GetSMstatus(isector)) {
      snprintf(cTagV,kTag,"BTRD%d",isector);
      TVirtualMC::GetMC()->Gspos("UTF1",1,cTagV,xpos, ypos,zpos,0,"ONLY");
      TVirtualMC::GetMC()->Gspos("UTF2",1,cTagV,xpos,-ypos,zpos,0,"ONLY");
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

  Int_t   ilayer = 0;

  Float_t xpos  = 0.0;
  Float_t ypos  = 0.0;
  Float_t zpos  = 0.0;

  const Int_t kTag = 100;
  Char_t  cTagV[kTag];
  Char_t  cTagM[kTag];

  const Int_t kNparTRD = 4;
  Float_t parTRD[kNparTRD];
  const Int_t kNparBOX = 3;
  Float_t parBOX[kNparBOX];
  const Int_t kNparTRP = 11;
  Float_t parTRP[kNparTRP];

  // The rotation matrices
  const Int_t kNmatrix = 7;
  Int_t   matrix[kNmatrix];
  TVirtualMC::GetMC()->Matrix(matrix[0], 100.0,   0.0,  90.0,  90.0,  10.0,   0.0);
  TVirtualMC::GetMC()->Matrix(matrix[1],  80.0,   0.0,  90.0,  90.0,  10.0, 180.0);
  TVirtualMC::GetMC()->Matrix(matrix[2],  90.0,   0.0,   0.0,   0.0,  90.0,  90.0);
  TVirtualMC::GetMC()->Matrix(matrix[3],  90.0, 180.0,   0.0, 180.0,  90.0,  90.0);
  TVirtualMC::GetMC()->Matrix(matrix[4], 170.0,   0.0,  80.0,   0.0,  90.0,  90.0);
  TVirtualMC::GetMC()->Matrix(matrix[5], 170.0, 180.0,  80.0, 180.0,  90.0,  90.0);
  TVirtualMC::GetMC()->Matrix(matrix[6], 180.0, 180.0,  90.0, 180.0,  90.0,  90.0);

  //
  // The carbon inserts in the top/bottom aluminum plates
  //

  const Int_t kNparCrb = 3;
  Float_t parCrb[kNparCrb];
  parCrb[0] = 0.0;
  parCrb[1] = 0.0;
  parCrb[2] = 0.0;
  TVirtualMC::GetMC()->Gsvolu("USCR","BOX ",idtmed[1326-1],parCrb,0);
  // Bottom 1 (all sectors)
  parCrb[0] =  77.49/2.0;
  parCrb[1] = 104.60/2.0;
  parCrb[2] = fgkSMpltT/2.0;
  xpos      =   0.0;
  ypos      =   0.0;
  zpos      = fgkSMpltT/2.0 - fgkSheight/2.0;
  TVirtualMC::GetMC()->Gsposp("USCR", 1,"UTS1", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR", 2,"UTS2", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR", 3,"UTS3", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR", 4,"UTS4", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  // Bottom 2 (all sectors)
  parCrb[0] =  77.49/2.0;
  parCrb[1] =  55.80/2.0;
  parCrb[2] = fgkSMpltT/2.0;
  xpos      =   0.0;
  ypos      =  85.6;
  zpos      = fgkSMpltT/2.0 - fgkSheight/2.0;
  TVirtualMC::GetMC()->Gsposp("USCR", 5,"UTS1", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR", 6,"UTS2", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR", 7,"UTS3", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR", 8,"UTS4", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR", 9,"UTS1", xpos,-ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",10,"UTS2", xpos,-ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",11,"UTS3", xpos,-ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",12,"UTS4", xpos,-ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  // Bottom 3 (all sectors)
  parCrb[0] =  77.49/2.0;
  parCrb[1] =  56.00/2.0;
  parCrb[2] = fgkSMpltT/2.0;
  xpos      =   0.0;
  ypos      = 148.5;
  zpos      = fgkSMpltT/2.0 - fgkSheight/2.0;
  TVirtualMC::GetMC()->Gsposp("USCR",13,"UTS1", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",14,"UTS2", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",15,"UTS3", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",16,"UTS4", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",17,"UTS1", xpos,-ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",18,"UTS2", xpos,-ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",19,"UTS3", xpos,-ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",20,"UTS4", xpos,-ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  // Bottom 4 (all sectors)
  parCrb[0] =  77.49/2.0;
  parCrb[1] = 118.00/2.0;
  parCrb[2] = fgkSMpltT/2.0;
  xpos      =   0.0;
  ypos      = 240.5;
  zpos      = fgkSMpltT/2.0 - fgkSheight/2.0;
  TVirtualMC::GetMC()->Gsposp("USCR",21,"UTS1", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",22,"UTS2", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",23,"UTS3", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",24,"UTS4", xpos, ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",25,"UTS1", xpos,-ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",26,"UTS2", xpos,-ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",27,"UTS3", xpos,-ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",28,"UTS4", xpos,-ypos, zpos,0,"ONLY",parCrb,kNparCrb);
  // Top 1 (only in front of PHOS)
  parCrb[0] = 111.48/2.0;
  parCrb[1] = 105.00/2.0;
  parCrb[2] = fgkSMpltT/2.0;
  xpos      =   0.0;
  ypos      =   0.0;
  zpos      = fgkSMpltT/2.0 - fgkSheight/2.0;
  TVirtualMC::GetMC()->Gsposp("USCR",29,"UTS2", xpos, ypos,-zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",30,"UTS3", xpos, ypos,-zpos,0,"ONLY",parCrb,kNparCrb);
  // Top 2 (only in front of PHOS)
  parCrb[0] = 111.48/2.0;
  parCrb[1] =  56.00/2.0;
  parCrb[2] = fgkSMpltT/2.0;
  xpos      =   0.0;
  ypos      =  85.5;
  zpos      = fgkSMpltT/2.0 - fgkSheight/2.0;
  TVirtualMC::GetMC()->Gsposp("USCR",31,"UTS2", xpos, ypos,-zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",32,"UTS3", xpos, ypos,-zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",33,"UTS2", xpos,-ypos,-zpos,0,"ONLY",parCrb,kNparCrb);
  TVirtualMC::GetMC()->Gsposp("USCR",34,"UTS3", xpos,-ypos,-zpos,0,"ONLY",parCrb,kNparCrb);

  //
  // The chamber support rails
  //

  const Float_t kSRLhgt  = 2.00;
  const Float_t kSRLwidA = 2.3;
  const Float_t kSRLwidB = 1.947;
  const Float_t kSRLdst  = 1.135;
  const Int_t   kNparSRL = 11;
  Float_t parSRL[kNparSRL];
  // Trapezoidal shape
  parSRL[ 0] = fgkSlength/2.0;
  parSRL[ 1] = 0.0;
  parSRL[ 2] = 0.0;
  parSRL[ 3] = kSRLhgt  /2.0;
  parSRL[ 4] = kSRLwidB /2.0;
  parSRL[ 5] = kSRLwidA /2.0;
  parSRL[ 6] = 5.0;
  parSRL[ 7] = kSRLhgt  /2.0;
  parSRL[ 8] = kSRLwidB /2.0;
  parSRL[ 9] = kSRLwidA /2.0;
  parSRL[10] = 5.0;
  TVirtualMC::GetMC()->Gsvolu("USRL","TRAP",idtmed[1301-1],parSRL,kNparSRL);

  xpos  = 0.0;
  ypos  = 0.0;
  zpos  = 0.0;
  for (ilayer = 1; ilayer < kNlayer; ilayer++) {
    xpos  = fgkCwidth[ilayer]/2.0 + kSRLwidA/2.0 + kSRLdst;
    ypos  = 0.0;
    zpos  = fgkVrocsm + fgkSMpltT - fgkCalZpos - fgkSheight/2.0  
          + fgkCraH + fgkCdrH - fgkCalH - kSRLhgt/2.0 
          + ilayer * (fgkCH + fgkVspace);
    TVirtualMC::GetMC()->Gspos("USRL",ilayer+1          ,"UTI1", xpos,ypos,zpos,matrix[2],"ONLY");
    TVirtualMC::GetMC()->Gspos("USRL",ilayer+1+  kNlayer,"UTI1",-xpos,ypos,zpos,matrix[3],"ONLY");
    TVirtualMC::GetMC()->Gspos("USRL",ilayer+1+2*kNlayer,"UTI2", xpos,ypos,zpos,matrix[2],"ONLY");
    TVirtualMC::GetMC()->Gspos("USRL",ilayer+1+3*kNlayer,"UTI2",-xpos,ypos,zpos,matrix[3],"ONLY");
    TVirtualMC::GetMC()->Gspos("USRL",ilayer+1+4*kNlayer,"UTI3", xpos,ypos,zpos,matrix[2],"ONLY");
    TVirtualMC::GetMC()->Gspos("USRL",ilayer+1+5*kNlayer,"UTI3",-xpos,ypos,zpos,matrix[3],"ONLY");
    TVirtualMC::GetMC()->Gspos("USRL",ilayer+1+6*kNlayer,"UTI4", xpos,ypos,zpos,matrix[2],"ONLY");
    TVirtualMC::GetMC()->Gspos("USRL",ilayer+1+7*kNlayer,"UTI4",-xpos,ypos,zpos,matrix[3],"ONLY");
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
  for (ilayer = 0; ilayer < kNlayer; ilayer++) {

    // The aluminum of the cross bars
    parSCB[0] = fgkCwidth[ilayer]/2.0 + kSRLdst/2.0;
    snprintf(cTagV,kTag,"USF%01d",ilayer);
    TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB);

    // The empty regions in the cross bars
    Float_t thkSCB = kSCBthk;
    if (ilayer < 2) {
      thkSCB *= 1.5;
    }
    parSCI[2] = parSCB[2] - thkSCB;
    parSCI[0] = parSCB[0]/4.0 - kSCBthk;
    snprintf(cTagV,kTag,"USI%01d",ilayer);
    TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parSCI,kNparSCI);

    snprintf(cTagV,kTag,"USI%01d",ilayer);
    snprintf(cTagM,kTag,"USF%01d",ilayer);
    ypos  = 0.0;
    zpos  = 0.0;
    xpos  =   parSCI[0] + thkSCB/2.0;
    TVirtualMC::GetMC()->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
    xpos  = - parSCI[0] - thkSCB/2.0;
    TVirtualMC::GetMC()->Gspos(cTagV,2,cTagM,xpos,ypos,zpos,0,"ONLY");
    xpos  =   3.0 * parSCI[0] + 1.5 * thkSCB;
    TVirtualMC::GetMC()->Gspos(cTagV,3,cTagM,xpos,ypos,zpos,0,"ONLY");
    xpos  = - 3.0 * parSCI[0] - 1.5 * thkSCB;
    TVirtualMC::GetMC()->Gspos(cTagV,4,cTagM,xpos,ypos,zpos,0,"ONLY");

    snprintf(cTagV,kTag,"USF%01d",ilayer);
    xpos  = 0.0;
    zpos  = fgkVrocsm + fgkSMpltT + parSCB[2] - fgkSheight/2.0 
          + ilayer * (fgkCH + fgkVspace);

    ypos  =   fgkClength[ilayer][2]/2.0 + fgkClength[ilayer][1];
    TVirtualMC::GetMC()->Gspos(cTagV, 1,"UTI1", xpos,ypos,zpos,0,"ONLY");
    TVirtualMC::GetMC()->Gspos(cTagV, 3,"UTI2", xpos,ypos,zpos,0,"ONLY");
    TVirtualMC::GetMC()->Gspos(cTagV, 5,"UTI3", xpos,ypos,zpos,0,"ONLY");
    TVirtualMC::GetMC()->Gspos(cTagV, 7,"UTI4", xpos,ypos,zpos,0,"ONLY");

    ypos  = - fgkClength[ilayer][2]/2.0 - fgkClength[ilayer][1];
    TVirtualMC::GetMC()->Gspos(cTagV, 2,"UTI1", xpos,ypos,zpos,0,"ONLY");
    TVirtualMC::GetMC()->Gspos(cTagV, 4,"UTI2", xpos,ypos,zpos,0,"ONLY");
    TVirtualMC::GetMC()->Gspos(cTagV, 6,"UTI3", xpos,ypos,zpos,0,"ONLY");
    TVirtualMC::GetMC()->Gspos(cTagV, 8,"UTI4", xpos,ypos,zpos,0,"ONLY");

  }

  //
  // The horizontal connections between the cross bars
  //

  const Int_t   kNparSCH = 3;
  Float_t parSCH[kNparSCH];

  for (ilayer = 1; ilayer < kNlayer-1; ilayer++) {

    parSCH[0] = fgkCwidth[ilayer]/2.0;
    parSCH[1] = (fgkClength[ilayer+1][2]/2.0 + fgkClength[ilayer+1][1]
               - fgkClength[ilayer  ][2]/2.0 - fgkClength[ilayer  ][1])/2.0;
    parSCH[2] = kSCHhgt/2.0;

    snprintf(cTagV,kTag,"USH%01d",ilayer);
    TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCH,kNparSCH);
    xpos  = 0.0;
    ypos  = fgkClength[ilayer][2]/2.0 + fgkClength[ilayer][1] + parSCH[1];
    zpos  = fgkVrocsm + fgkSMpltT - kSCHhgt/2.0 - fgkSheight/2.0 
          + (ilayer+1) * (fgkCH + fgkVspace);
    TVirtualMC::GetMC()->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY");
    TVirtualMC::GetMC()->Gspos(cTagV,3,"UTI2", xpos,ypos,zpos,0,"ONLY");
    TVirtualMC::GetMC()->Gspos(cTagV,5,"UTI3", xpos,ypos,zpos,0,"ONLY");
    TVirtualMC::GetMC()->Gspos(cTagV,7,"UTI4", xpos,ypos,zpos,0,"ONLY");
    ypos  = -ypos;
    TVirtualMC::GetMC()->Gspos(cTagV,2,"UTI1", xpos,ypos,zpos,0,"ONLY");
    TVirtualMC::GetMC()->Gspos(cTagV,4,"UTI2", xpos,ypos,zpos,0,"ONLY");
    TVirtualMC::GetMC()->Gspos(cTagV,6,"UTI3", xpos,ypos,zpos,0,"ONLY");
    TVirtualMC::GetMC()->Gspos(cTagV,8,"UTI4", xpos,ypos,zpos,0,"ONLY");

  }

  //
  // The aymmetric flat frame in the middle
  //

  // The envelope volume (aluminum)
  parTRD[0]  =  87.60/2.0;
  parTRD[1]  = 114.00/2.0;
  parTRD[2]  =   1.20/2.0;
  parTRD[3]  =  71.30/2.0;
  TVirtualMC::GetMC()->Gsvolu("USDB","TRD1",idtmed[1301-1],parTRD,kNparTRD);
  // Empty spaces (air)
  parTRP[ 0] =   1.20/2.0;
  parTRP[ 1] =   0.0;
  parTRP[ 2] =   0.0;
  parTRP[ 3] =  27.00/2.0;
  parTRP[ 4] =  50.60/2.0;
  parTRP[ 5] =   5.00/2.0;
  parTRP[ 6] =   3.5;
  parTRP[ 7] =  27.00/2.0;
  parTRP[ 8] =  50.60/2.0;
  parTRP[ 9] =   5.00/2.0;
  parTRP[10] =   3.5;
  TVirtualMC::GetMC()->Gsvolu("USD1","TRAP",idtmed[1302-1],parTRP,kNparTRP);
  xpos       =  18.0;
  ypos       =   0.0;
  zpos       =   27.00/2.0 - 71.3/2.0;
  TVirtualMC::GetMC()->Gspos("USD1",1,"USDB", xpos, ypos, zpos,matrix[2],"ONLY");
  // Empty spaces (air)
  parTRP[ 0] =   1.20/2.0;
  parTRP[ 1] =   0.0;
  parTRP[ 2] =   0.0;
  parTRP[ 3] =  33.00/2.0;
  parTRP[ 4] =   5.00/2.0;
  parTRP[ 5] =  62.10/2.0;
  parTRP[ 6] =   3.5;
  parTRP[ 7] =  33.00/2.0;
  parTRP[ 8] =   5.00/2.0;
  parTRP[ 9] =  62.10/2.0;
  parTRP[10] =   3.5;
  TVirtualMC::GetMC()->Gsvolu("USD2","TRAP",idtmed[1302-1],parTRP,kNparTRP);
  xpos       =  21.0;
  ypos       =   0.0;
  zpos       =  71.3/2.0 - 33.0/2.0;
  TVirtualMC::GetMC()->Gspos("USD2",1,"USDB", xpos, ypos, zpos,matrix[2],"ONLY");
  // Empty spaces (air)
  parBOX[ 0] =  22.50/2.0;
  parBOX[ 1] =   1.20/2.0;
  parBOX[ 2] =  70.50/2.0;
  TVirtualMC::GetMC()->Gsvolu("USD3","BOX ",idtmed[1302-1],parBOX,kNparBOX);
  xpos       = -25.75;
  ypos       =   0.0;
  zpos       =   0.4;
  TVirtualMC::GetMC()->Gspos("USD3",1,"USDB", xpos, ypos, zpos,        0,"ONLY");
  // Empty spaces (air)
  parTRP[ 0] =   1.20/2.0;
  parTRP[ 1] =   0.0;
  parTRP[ 2] =   0.0;
  parTRP[ 3] =  25.50/2.0;
  parTRP[ 4] =   5.00/2.0;
  parTRP[ 5] =  65.00/2.0;
  parTRP[ 6] =  -1.0;
  parTRP[ 7] =  25.50/2.0;
  parTRP[ 8] =   5.00/2.0;
  parTRP[ 9] =  65.00/2.0;
  parTRP[10] =  -1.0;
  TVirtualMC::GetMC()->Gsvolu("USD4","TRAP",idtmed[1302-1],parTRP,kNparTRP);
  xpos       =   2.0;
  ypos       =   0.0;
  zpos       =  -1.6;
  TVirtualMC::GetMC()->Gspos("USD4",1,"USDB", xpos, ypos, zpos,matrix[6],"ONLY");
  // Empty spaces (air)
  parTRP[ 0] =   1.20/2.0;
  parTRP[ 1] =   0.0;
  parTRP[ 2] =   0.0;
  parTRP[ 3] =  23.50/2.0;
  parTRP[ 4] =  63.50/2.0;
  parTRP[ 5] =   5.00/2.0;
  parTRP[ 6] =  16.0;
  parTRP[ 7] =  23.50/2.0;
  parTRP[ 8] =  63.50/2.0;
  parTRP[ 9] =   5.00/2.0;
  parTRP[10] =  16.0;
  TVirtualMC::GetMC()->Gsvolu("USD5","TRAP",idtmed[1302-1],parTRP,kNparTRP);
  xpos       =  36.5;
  ypos       =   0.0;
  zpos       =  -1.5;
  TVirtualMC::GetMC()->Gspos("USD5",1,"USDB", xpos, ypos, zpos,matrix[5],"ONLY");
  // Empty spaces (air)
  parTRP[ 0] =   1.20/2.0;
  parTRP[ 1] =   0.0;
  parTRP[ 2] =   0.0;
  parTRP[ 3] =  70.50/2.0;
  parTRP[ 4] =   4.50/2.0;
  parTRP[ 5] =  16.50/2.0;
  parTRP[ 6] =  -5.0;
  parTRP[ 7] =  70.50/2.0;
  parTRP[ 8] =   4.50/2.0;
  parTRP[ 9] =  16.50/2.0;
  parTRP[10] =  -5.0;
  TVirtualMC::GetMC()->Gsvolu("USD6","TRAP",idtmed[1302-1],parTRP,kNparTRP);
  xpos       = -43.7;
  ypos       =   0.0;
  zpos       =   0.4;
  TVirtualMC::GetMC()->Gspos("USD6",1,"USDB", xpos, ypos, zpos,matrix[2],"ONLY");
  xpos       =   0.0;
  ypos       =   fgkClength[5][2]/2.0;
  zpos       =   0.04;
  TVirtualMC::GetMC()->Gspos("USDB",1,"UTI1", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USDB",2,"UTI1", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USDB",3,"UTI2", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USDB",4,"UTI2", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USDB",5,"UTI3", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USDB",6,"UTI3", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USDB",7,"UTI4", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USDB",8,"UTI4", xpos,-ypos, zpos,        0,"ONLY");
  // Upper bar (aluminum)
  parBOX[0] = 95.00/2.0;
  parBOX[1] =  1.20/2.0;
  parBOX[2] =  3.00/2.0;
  TVirtualMC::GetMC()->Gsvolu("USD7","BOX ",idtmed[1301-1],parBOX,kNparBOX);
  xpos       =   0.0;
  ypos       =   fgkClength[5][2]/2.0;
  zpos       =   fgkSheight/2.0 - fgkSMpltT  - 3.00/2.0;
  TVirtualMC::GetMC()->Gspos("USD7",1,"UTI1", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD7",2,"UTI1", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD7",3,"UTI2", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD7",4,"UTI2", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD7",5,"UTI3", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD7",6,"UTI3", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD7",7,"UTI4", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD7",8,"UTI4", xpos,-ypos, zpos,        0,"ONLY");
  // Lower bar (aluminum)
  parBOX[0] = 90.22/2.0;
  parBOX[1] =  1.20/2.0;
  parBOX[2] =  1.74/2.0;
  TVirtualMC::GetMC()->Gsvolu("USD8","BOX ",idtmed[1301-1],parBOX,kNparBOX);
  xpos       =   0.0;
  ypos       =   fgkClength[5][2]/2.0 - 0.1;
  zpos       =  -fgkSheight/2.0 + fgkSMpltT + 2.27;
  TVirtualMC::GetMC()->Gspos("USD8",1,"UTI1", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD8",2,"UTI1", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD8",3,"UTI2", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD8",4,"UTI2", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD8",5,"UTI3", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD8",6,"UTI3", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD8",7,"UTI4", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD8",8,"UTI4", xpos,-ypos, zpos,        0,"ONLY");
  // Lower bar (aluminum)
  parBOX[0] = 82.60/2.0;
  parBOX[1] =  1.20/2.0;
  parBOX[2] =  1.40/2.0;
  TVirtualMC::GetMC()->Gsvolu("USD9","BOX ",idtmed[1301-1],parBOX,kNparBOX);
  xpos       =   0.0;
  ypos       =   fgkClength[5][2]/2.0;
  zpos       =  -fgkSheight/2.0 + fgkSMpltT + 1.40/2.0;
  TVirtualMC::GetMC()->Gspos("USD9",1,"UTI1", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD9",2,"UTI1", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD9",3,"UTI2", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD9",4,"UTI2", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD9",5,"UTI3", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD9",6,"UTI3", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD9",7,"UTI4", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USD9",8,"UTI4", xpos,-ypos, zpos,        0,"ONLY");
  // Front sheet (aluminum)
  parTRP[ 0] =   0.10/2.0;
  parTRP[ 1] =   0.0;
  parTRP[ 2] =   0.0;
  parTRP[ 3] =  74.50/2.0;
  parTRP[ 4] =  31.70/2.0;
  parTRP[ 5] =  44.00/2.0;
  parTRP[ 6] =  -5.0;
  parTRP[ 7] =  74.50/2.0;
  parTRP[ 8] =  31.70/2.0;
  parTRP[ 9] =  44.00/2.0;
  parTRP[10] =  -5.0;
  TVirtualMC::GetMC()->Gsvolu("USDF","TRAP",idtmed[1302-1],parTRP,kNparTRP);
  xpos       = -32.0;
  ypos       =   fgkClength[5][2]/2.0 + 1.20/2.0 + 0.10/2.0;
  zpos       =   0.0;
  TVirtualMC::GetMC()->Gspos("USDF",1,"UTI1", xpos, ypos, zpos,matrix[2],"ONLY");
  TVirtualMC::GetMC()->Gspos("USDF",2,"UTI1", xpos,-ypos, zpos,matrix[2],"ONLY");
  TVirtualMC::GetMC()->Gspos("USDF",3,"UTI2", xpos, ypos, zpos,matrix[2],"ONLY");
  TVirtualMC::GetMC()->Gspos("USDF",4,"UTI2", xpos,-ypos, zpos,matrix[2],"ONLY");
  TVirtualMC::GetMC()->Gspos("USDF",5,"UTI3", xpos, ypos, zpos,matrix[2],"ONLY");
  TVirtualMC::GetMC()->Gspos("USDF",6,"UTI3", xpos,-ypos, zpos,matrix[2],"ONLY");
  TVirtualMC::GetMC()->Gspos("USDF",7,"UTI4", xpos, ypos, zpos,matrix[2],"ONLY");
  TVirtualMC::GetMC()->Gspos("USDF",8,"UTI4", xpos,-ypos, zpos,matrix[2],"ONLY");

  //
  // The flat frame in front of the chambers
  //

  // The envelope volume (aluminum)
  parTRD[0]  =  90.00/2.0 - 0.1;
  parTRD[1]  = 114.00/2.0 - 0.1;
  parTRD[2]  =   1.50/2.0;
  parTRD[3]  =  70.30/2.0;
  TVirtualMC::GetMC()->Gsvolu("USCB","TRD1",idtmed[1301-1],parTRD,kNparTRD);
  // Empty spaces (air)
  parTRD[0]  =  87.00/2.0;
  parTRD[1]  =  10.00/2.0;
  parTRD[2]  =   1.50/2.0;
  parTRD[3]  =  26.35/2.0;
  TVirtualMC::GetMC()->Gsvolu("USC1","TRD1",idtmed[1302-1],parTRD,kNparTRD);
  xpos       =  0.0;
  ypos       =  0.0;
  zpos       = 26.35/2.0 - 70.3/2.0;
  TVirtualMC::GetMC()->Gspos("USC1",1,"USCB",xpos,ypos,zpos,0,"ONLY");
  // Empty spaces (air)
  parTRD[0]  =  10.00/2.0;
  parTRD[1]  = 111.00/2.0;
  parTRD[2]  =   1.50/2.0;
  parTRD[3]  =  35.05/2.0;
  TVirtualMC::GetMC()->Gsvolu("USC2","TRD1",idtmed[1302-1],parTRD,kNparTRD);
  xpos       =  0.0;
  ypos       =  0.0;
  zpos       = 70.3/2.0 - 35.05/2.0;
  TVirtualMC::GetMC()->Gspos("USC2",1,"USCB",xpos,ypos,zpos,0,"ONLY");
  // Empty spaces (air)
  parTRP[ 0] =   1.50/2.0;
  parTRP[ 1] =   0.0;
  parTRP[ 2] =   0.0;
  parTRP[ 3] =  37.60/2.0;
  parTRP[ 4] =  63.90/2.0;
  parTRP[ 5] =   8.86/2.0;
  parTRP[ 6] =  16.0;
  parTRP[ 7] =  37.60/2.0;
  parTRP[ 8] =  63.90/2.0;
  parTRP[ 9] =   8.86/2.0;
  parTRP[10] =  16.0;
  TVirtualMC::GetMC()->Gsvolu("USC3","TRAP",idtmed[1302-1],parTRP,kNparTRP);
  xpos       = -30.5;
  ypos       =   0.0;
  zpos       =  -2.0;
  TVirtualMC::GetMC()->Gspos("USC3",1,"USCB", xpos, ypos, zpos,matrix[4],"ONLY");
  TVirtualMC::GetMC()->Gspos("USC3",2,"USCB",-xpos, ypos, zpos,matrix[5],"ONLY");
  xpos       =   0.0;
  ypos       =   fgkClength[5][2]/2.0 + fgkClength[5][1] + fgkClength[5][0];
  zpos       =   0.0;
  TVirtualMC::GetMC()->Gspos("USCB",1,"UTI1", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USCB",2,"UTI1", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USCB",3,"UTI2", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USCB",4,"UTI2", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USCB",5,"UTI3", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USCB",6,"UTI3", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USCB",7,"UTI4", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USCB",8,"UTI4", xpos,-ypos, zpos,        0,"ONLY");
  // Upper bar (aluminum)
  parBOX[0] = 95.00/2.0;
  parBOX[1] =  1.50/2.0;
  parBOX[2] =  3.00/2.0;
  TVirtualMC::GetMC()->Gsvolu("USC4","BOX ",idtmed[1301-1],parBOX,kNparBOX);
  xpos       =   0.0;
  ypos       =   fgkClength[5][2]/2.0 + fgkClength[5][1] + fgkClength[5][0];
  zpos       =   fgkSheight/2.0 - fgkSMpltT - 3.00/2.0;
  TVirtualMC::GetMC()->Gspos("USC4",1,"UTI1", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC4",2,"UTI1", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC4",3,"UTI2", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC4",4,"UTI2", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC4",5,"UTI3", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC4",6,"UTI3", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC4",7,"UTI4", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC4",8,"UTI4", xpos,-ypos, zpos,        0,"ONLY");
  // Lower bar (aluminum)
  parBOX[0] = 90.22/2.0;
  parBOX[1] =  1.50/2.0;
  parBOX[2] =  2.00/2.0;
  TVirtualMC::GetMC()->Gsvolu("USC5","BOX ",idtmed[1301-1],parBOX,kNparBOX);
  xpos       =   0.0;
  ypos       =   fgkClength[5][2]/2.0 + fgkClength[5][1] + fgkClength[5][0];
  zpos       =  -fgkSheight/2.0 + fgkSMpltT + 2.60;
  TVirtualMC::GetMC()->Gspos("USC5",1,"UTI1", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC5",2,"UTI1", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC5",3,"UTI2", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC5",4,"UTI2", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC5",5,"UTI3", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC5",6,"UTI3", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC5",7,"UTI4", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC5",8,"UTI4", xpos,-ypos, zpos,        0,"ONLY");
  // Lower bar (aluminum)
  parBOX[0] = 82.60/2.0;
  parBOX[1] =  1.50/2.0;
  parBOX[2] =  1.60/2.0;
  TVirtualMC::GetMC()->Gsvolu("USC6","BOX ",idtmed[1301-1],parBOX,kNparBOX);
  xpos       =   0.0;
  ypos       =   fgkClength[5][2]/2.0 + fgkClength[5][1] + fgkClength[5][0];
  zpos       =  -fgkSheight/2.0 + fgkSMpltT + 1.60/2.0;
  TVirtualMC::GetMC()->Gspos("USC6",1,"UTI1", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC6",2,"UTI1", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC6",3,"UTI2", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC6",4,"UTI2", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC6",5,"UTI3", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC6",6,"UTI3", xpos,-ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC6",7,"UTI4", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USC6",8,"UTI4", xpos,-ypos, zpos,        0,"ONLY");

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
  const Float_t kSCLposzUa = 1.65;
  const Float_t kSCLposzUb = 0.3;
  // Vertical
  parSCL[0]  = kSCLthkUa /2.0;
  parSCL[1]  = fgkSlength/2.0;
  parSCL[2]  = kSCLwidUa /2.0;
  TVirtualMC::GetMC()->Gsvolu("USL1","BOX ",idtmed[1301-1],parSCL,kNparSCL);
  xpos  =   fgkSwidth2/2.0 - fgkSMpltT - kSCLposxUa;
  ypos  =   0.0;
  zpos  =   fgkSheight/2.0 - fgkSMpltT - kSCLposzUa;
  TVirtualMC::GetMC()->Gspos("USL1",1,"UTI1", xpos,ypos,zpos,matrix[0],"ONLY");
  TVirtualMC::GetMC()->Gspos("USL1",3,"UTI4", xpos,ypos,zpos,matrix[0],"ONLY");
  xpos  = -xpos;
  TVirtualMC::GetMC()->Gspos("USL1",2,"UTI1", xpos,ypos,zpos,matrix[1],"ONLY");
  TVirtualMC::GetMC()->Gspos("USL1",4,"UTI4", xpos,ypos,zpos,matrix[1],"ONLY");
  // Horizontal
  parSCL[0]  = kSCLwidUb /2.0;
  parSCL[1]  = fgkSlength/2.0;
  parSCL[2]  = kSCLthkUb /2.0;
  TVirtualMC::GetMC()->Gsvolu("USL2","BOX ",idtmed[1301-1],parSCL,kNparSCL);
  xpos  =   fgkSwidth2/2.0 - fgkSMpltT - kSCLposxUb;
  ypos  =   0.0;
  zpos  =   fgkSheight/2.0 - fgkSMpltT - kSCLposzUb; 
  TVirtualMC::GetMC()->Gspos("USL2",1,"UTI1", xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USL2",3,"UTI2", xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USL2",5,"UTI3", xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USL2",7,"UTI4", xpos,ypos,zpos,        0,"ONLY");
  xpos  = -xpos;
  TVirtualMC::GetMC()->Gspos("USL2",2,"UTI1", xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USL2",4,"UTI2", xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USL2",6,"UTI3", xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USL2",8,"UTI4", xpos,ypos,zpos,        0,"ONLY");

  // Lower ledges 
  // Thickness of the corner ledges
  const Float_t kSCLthkLa  =  2.464; 
  const Float_t kSCLthkLb  =  1.0; 
  // Width of the corner ledges
  const Float_t kSCLwidLa  =  8.3;
  const Float_t kSCLwidLb  =  4.0;
  // Position of the corner ledges
  const Float_t kSCLposxLa = (3.0 * kSCLthkLb - kSCLthkLa) / 4.0 + 0.05;
  const Float_t kSCLposxLb = kSCLthkLb + kSCLwidLb/2.0 + 0.05;
  const Float_t kSCLposzLa = kSCLwidLa/2.0;
  const Float_t kSCLposzLb = kSCLthkLb/2.0;
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
  TVirtualMC::GetMC()->Gsvolu("USL3","TRAP",idtmed[1301-1],parSCLb,kNparSCLb);
  xpos  =   fgkSwidth1/2.0 - fgkSMpltT - kSCLposxLa;
  ypos  =   0.0;
  zpos  = - fgkSheight/2.0 + fgkSMpltT + kSCLposzLa;
  TVirtualMC::GetMC()->Gspos("USL3",1,"UTI1", xpos,ypos,zpos,matrix[2],"ONLY");
  TVirtualMC::GetMC()->Gspos("USL3",3,"UTI2", xpos,ypos,zpos,matrix[2],"ONLY");
  TVirtualMC::GetMC()->Gspos("USL3",5,"UTI3", xpos,ypos,zpos,matrix[2],"ONLY");
  TVirtualMC::GetMC()->Gspos("USL3",7,"UTI4", xpos,ypos,zpos,matrix[2],"ONLY");
  xpos  = -xpos;
  TVirtualMC::GetMC()->Gspos("USL3",2,"UTI1", xpos,ypos,zpos,matrix[3],"ONLY");
  TVirtualMC::GetMC()->Gspos("USL3",4,"UTI2", xpos,ypos,zpos,matrix[3],"ONLY");
  TVirtualMC::GetMC()->Gspos("USL3",6,"UTI3", xpos,ypos,zpos,matrix[3],"ONLY");
  TVirtualMC::GetMC()->Gspos("USL3",8,"UTI4", xpos,ypos,zpos,matrix[3],"ONLY");
  // Horizontal part
  parSCL[0]  = kSCLwidLb /2.0;
  parSCL[1]  = fgkSlength/2.0;
  parSCL[2]  = kSCLthkLb /2.0;
  TVirtualMC::GetMC()->Gsvolu("USL4","BOX ",idtmed[1301-1],parSCL,kNparSCL);
  xpos  =   fgkSwidth1/2.0 - fgkSMpltT - kSCLposxLb;
  ypos  =   0.0;
  zpos  = - fgkSheight/2.0 + fgkSMpltT + kSCLposzLb;
  TVirtualMC::GetMC()->Gspos("USL4",1,"UTI1", xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USL4",3,"UTI2", xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USL4",5,"UTI3", xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USL4",7,"UTI4", xpos,ypos,zpos,        0,"ONLY");
  xpos  = -xpos;
  TVirtualMC::GetMC()->Gspos("USL4",2,"UTI1", xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USL4",4,"UTI2", xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USL4",6,"UTI3", xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("USL4",8,"UTI4", xpos,ypos,zpos,        0,"ONLY");

  //
  // Aluminum plates in the front part of the super modules
  //

  const Int_t kNparTrd = 4;
  Float_t parTrd[kNparTrd];
  parTrd[0] = fgkSwidth1/2.0 - 2.5;
  parTrd[1] = fgkSwidth2/2.0 - 2.5;
  parTrd[2] = fgkSMpltT /2.0;
  parTrd[3] = fgkSheight/2.0 - 1.0;
  TVirtualMC::GetMC()->Gsvolu("UTA1","TRD1",idtmed[1301-1],parTrd,kNparTrd);
  xpos      =  0.0;
  ypos      =  fgkSMpltT/2.0 - fgkFlength/2.0;
  zpos      = -0.5;
  TVirtualMC::GetMC()->Gspos("UTA1",1,"UTF1",xpos, ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTA1",2,"UTF2",xpos,-ypos,zpos,        0,"ONLY");

  const Int_t kNparPlt = 3;
  Float_t parPlt[kNparPlt];
  parPlt[0] =  0.0;
  parPlt[1] =  0.0;
  parPlt[2] =  0.0;
  TVirtualMC::GetMC()->Gsvolu("UTA2","BOX ",idtmed[1301-1],parPlt,0);
  xpos      =  0.0;
  ypos      =  0.0;
  zpos      =  fgkSheight/2.0 - fgkSMpltT/2.0;
  parPlt[0] = fgkSwidth2/2.0 - 0.2;
  parPlt[1] = fgkFlength/2.0;
  parPlt[2] = fgkSMpltT /2.0;
  TVirtualMC::GetMC()->Gsposp("UTA2",1,"UTF2",xpos,ypos,zpos
                    ,        0,"ONLY",parPlt,kNparPlt);
  xpos      = (fgkSwidth1 + fgkSwidth2)/4.0 - fgkSMpltT/2.0 - 0.0016;
  ypos      =  0.0;
  zpos      =  0.0;
  parPlt[0] = fgkSMpltT /2.0;
  parPlt[1] = fgkFlength/2.0;
  parPlt[2] = fgkSheight/2.0;
  TVirtualMC::GetMC()->Gsposp("UTA2",2,"UTF2", xpos,ypos,zpos
                    ,matrix[0],"ONLY",parPlt,kNparPlt);
  TVirtualMC::GetMC()->Gsposp("UTA2",3,"UTF2",-xpos,ypos,zpos
                    ,matrix[1],"ONLY",parPlt,kNparPlt);

  // Additional aluminum bar
  parBOX[0] = 80.0/2.0;
  parBOX[1] =  1.0/2.0;
  parBOX[2] = 10.0/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTA3","BOX ",idtmed[1301-1],parBOX,kNparBOX);
  xpos      =  0.0;
  ypos      =  1.0/2.0 + fgkSMpltT - fgkFlength/2.0;
  zpos      =  fgkSheight/2.0 - 1.5 - 10.0/2.0;
  TVirtualMC::GetMC()->Gspos("UTA3",1,"UTF1", xpos, ypos, zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTA3",2,"UTF2", xpos,-ypos, zpos,        0,"ONLY");

}

//_____________________________________________________________________________
void AliTRDgeometry::CreateServices(Int_t *idtmed)
{
  //
  // Create the geometry of the services
  //
  // Names of the TRD services volumina
  //
  //        UTC1    Cooling arterias (Al)
  //        UTC2    Cooling arterias (Water)
  //        UUxx    Volumes for the services at the chambers (Air)
  //        UMCM    Readout MCMs     (G10/Cu/Si)
  //        UDCS    DCSs boards      (G10/Cu)
  //        UTP1    Power bars       (Cu)
  //        UTCP    Cooling pipes    (Fe)
  //        UTCH    Cooling pipes    (Water)
  //        UTPL    Power lines      (Cu)
  //        UTGD    Gas distribution box (V2A)
  //

  Int_t   ilayer = 0;
  Int_t   istack = 0;

  Float_t xpos  = 0.0;
  Float_t ypos  = 0.0;
  Float_t zpos  = 0.0;

  const Int_t kTag = 100;
  Char_t  cTagV[kTag];

  const Int_t kNparBox  = 3;
  Float_t parBox[kNparBox];

  const Int_t kNparTube = 3;
  Float_t parTube[kNparTube];

  // Services inside the baby frame
  const Float_t kBBMdz = 223.0;
  const Float_t kBBSdz =   8.5;

  // Services inside the back frame
  const Float_t kBFMdz = 118.0;
  const Float_t kBFSdz =   8.5;

  // The rotation matrices
  const Int_t kNmatrix = 10;
  Int_t   matrix[kNmatrix];
  TVirtualMC::GetMC()->Matrix(matrix[0], 100.0,   0.0,  90.0,  90.0,  10.0,   0.0); // rotation around y-axis
  TVirtualMC::GetMC()->Matrix(matrix[1],  80.0,   0.0,  90.0,  90.0,  10.0, 180.0); // rotation around y-axis
  TVirtualMC::GetMC()->Matrix(matrix[2],   0.0,   0.0,  90.0,  90.0,  90.0,   0.0);
  TVirtualMC::GetMC()->Matrix(matrix[3], 180.0,   0.0,  90.0,  90.0,  90.0, 180.0);
  TVirtualMC::GetMC()->Matrix(matrix[4],  90.0,   0.0,   0.0,   0.0,  90.0,  90.0);
  TVirtualMC::GetMC()->Matrix(matrix[5], 100.0,   0.0,  90.0, 270.0,  10.0,   0.0);
  TVirtualMC::GetMC()->Matrix(matrix[6],  80.0,   0.0,  90.0, 270.0,  10.0, 180.0);
  TVirtualMC::GetMC()->Matrix(matrix[7],  90.0,  10.0,  90.0, 100.0,   0.0,   0.0); // rotation around z-axis
  TVirtualMC::GetMC()->Matrix(matrix[8],  90.0, 350.0,  90.0,  80.0,   0.0,   0.0); // rotation around z-axis
  TVirtualMC::GetMC()->Matrix(matrix[9],  90.0,  90.0,  90.0, 180.0,   0.0,   0.0); // rotation around z-axis
    
  //
  // The cooling arterias
  //

  // Width of the cooling arterias
  const Float_t kCOLwid  =  0.8; 
  // Height of the cooling arterias
  const Float_t kCOLhgt  =  6.5;
  // Positioning of the cooling 
  const Float_t kCOLposx =  1.0;
  const Float_t kCOLposz = -1.2;
  // Thickness of the walls of the cooling arterias
  const Float_t kCOLthk  =  0.1;
  const Int_t   kNparCOL =  3;
  Float_t parCOL[kNparCOL];
  parCOL[0] = 0.0;
  parCOL[1] = 0.0;
  parCOL[2] = 0.0;
  TVirtualMC::GetMC()->Gsvolu("UTC1","BOX ",idtmed[1308-1],parCOL,0);
  TVirtualMC::GetMC()->Gsvolu("UTC3","BOX ",idtmed[1308-1],parCOL,0);
  parCOL[0] =  kCOLwid/2.0 - kCOLthk;
  parCOL[1] = -1.0;
  parCOL[2] =  kCOLhgt/2.0 - kCOLthk;
  TVirtualMC::GetMC()->Gsvolu("UTC2","BOX ",idtmed[1314-1],parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsvolu("UTC4","BOX ",idtmed[1314-1],parCOL,kNparCOL);

  xpos  = 0.0;
  ypos  = 0.0;
  zpos  = 0.0;
  TVirtualMC::GetMC()->Gspos("UTC2",1,"UTC1", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTC4",1,"UTC3", xpos,ypos,zpos,0,"ONLY");

  for (ilayer = 1; ilayer < kNlayer; ilayer++) { 

    // Along the chambers
    xpos      = fgkCwidth[ilayer]/2.0 + kCOLwid/2.0 + kCOLposx;
    ypos      = 0.0;
    zpos      = fgkVrocsm + fgkSMpltT - fgkCalZpos 
              + kCOLhgt/2.0 - fgkSheight/2.0 + kCOLposz 
              + ilayer * (fgkCH + fgkVspace);
    parCOL[0] = kCOLwid   /2.0;
    parCOL[1] = fgkSlength/2.0;
    parCOL[2] = kCOLhgt   /2.0;
    TVirtualMC::GetMC()->Gsposp("UTC1",ilayer           ,"UTI1", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parCOL,kNparCOL);
    TVirtualMC::GetMC()->Gsposp("UTC1",ilayer+   kNlayer,"UTI1",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parCOL,kNparCOL);
    TVirtualMC::GetMC()->Gsposp("UTC1",ilayer+ 6*kNlayer,"UTI2", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parCOL,kNparCOL);
    TVirtualMC::GetMC()->Gsposp("UTC1",ilayer+ 7*kNlayer,"UTI2",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parCOL,kNparCOL);
    TVirtualMC::GetMC()->Gsposp("UTC1",ilayer+ 8*kNlayer,"UTI3", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parCOL,kNparCOL);
    TVirtualMC::GetMC()->Gsposp("UTC1",ilayer+ 9*kNlayer,"UTI3",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parCOL,kNparCOL);
    TVirtualMC::GetMC()->Gsposp("UTC1",ilayer+10*kNlayer,"UTI4", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parCOL,kNparCOL);
    TVirtualMC::GetMC()->Gsposp("UTC1",ilayer+11*kNlayer,"UTI4",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parCOL,kNparCOL);

    // Front of supermodules
    xpos      = fgkCwidth[ilayer]/2.0 + kCOLwid/2.0 + kCOLposx;
    ypos      = 0.0;
    zpos      = fgkVrocsm + fgkSMpltT - fgkCalZpos 
              + kCOLhgt/2.0 - fgkSheight/2.0 + kCOLposz 
              + ilayer * (fgkCH + fgkVspace);
    parCOL[0] = kCOLwid   /2.0;
    parCOL[1] = fgkFlength/2.0;
    parCOL[2] = kCOLhgt   /2.0;
    TVirtualMC::GetMC()->Gsposp("UTC3",ilayer+2*kNlayer,"UTF1", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parCOL,kNparCOL);
    TVirtualMC::GetMC()->Gsposp("UTC3",ilayer+3*kNlayer,"UTF1",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parCOL,kNparCOL);
    TVirtualMC::GetMC()->Gsposp("UTC3",ilayer+4*kNlayer,"UTF2", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parCOL,kNparCOL);
    TVirtualMC::GetMC()->Gsposp("UTC3",ilayer+5*kNlayer,"UTF2",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parCOL,kNparCOL);

  }

  for (ilayer = 1; ilayer < kNlayer; ilayer++) { 

    // In baby frame
    xpos      = fgkCwidth[ilayer]/2.0 + kCOLwid/2.0 + kCOLposx - 2.5;
    ypos      = kBBSdz/2.0 - kBBMdz/2.0;
    zpos      = fgkVrocsm + fgkSMpltT - fgkCalZpos 
              + kCOLhgt/2.0 - fgkSheight/2.0 + kCOLposz 
              + ilayer * (fgkCH + fgkVspace);
    parCOL[0] = kCOLwid/2.0;
    parCOL[1] = kBBSdz /2.0;
    parCOL[2] = kCOLhgt/2.0;
    TVirtualMC::GetMC()->Gsposp("UTC3",ilayer+6*kNlayer,"BBTRD", xpos, ypos, zpos
                      ,matrix[0],"ONLY",parCOL,kNparCOL);
    TVirtualMC::GetMC()->Gsposp("UTC3",ilayer+7*kNlayer,"BBTRD",-xpos, ypos, zpos
                      ,matrix[1],"ONLY",parCOL,kNparCOL);

  }

  for (ilayer = 1; ilayer < kNlayer; ilayer++) { 

    // In back frame
    xpos      = fgkCwidth[ilayer]/2.0 + kCOLwid/2.0 + kCOLposx - 0.3;
    ypos      = -kBFSdz/2.0 + kBFMdz/2.0;
    zpos      = fgkVrocsm + fgkSMpltT - fgkCalZpos 
              + kCOLhgt/2.0 - fgkSheight/2.0 + kCOLposz 
              + ilayer * (fgkCH + fgkVspace);
    parCOL[0] = kCOLwid/2.0;
    parCOL[1] = kBFSdz /2.0;
    parCOL[2] = kCOLhgt/2.0;
    TVirtualMC::GetMC()->Gsposp("UTC3",ilayer+6*kNlayer,"BFTRD", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parCOL,kNparCOL);
    TVirtualMC::GetMC()->Gsposp("UTC3",ilayer+7*kNlayer,"BFTRD",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parCOL,kNparCOL);

  }

  // The upper most layer
  // Along the chambers
  xpos      = fgkCwidth[5]/2.0 - kCOLhgt/2.0 - 1.3;
  ypos      = 0.0;
  zpos      = fgkSheight/2.0 - fgkSMpltT - 0.4 - kCOLwid/2.0; 
  parCOL[0] = kCOLwid   /2.0;
  parCOL[1] = fgkSlength/2.0;
  parCOL[2] = kCOLhgt   /2.0;
  TVirtualMC::GetMC()->Gsposp("UTC1",6           ,"UTI1", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsposp("UTC1",6+   kNlayer,"UTI1",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsposp("UTC1",6+ 6*kNlayer,"UTI2", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsposp("UTC1",6+ 7*kNlayer,"UTI2",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsposp("UTC1",6+ 8*kNlayer,"UTI3", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsposp("UTC1",6+ 9*kNlayer,"UTI3",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsposp("UTC1",6+10*kNlayer,"UTI4", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsposp("UTC1",6+11*kNlayer,"UTI4",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  // Front of supermodules
  xpos      = fgkCwidth[5]/2.0 - kCOLhgt/2.0 - 1.3;
  ypos      = 0.0;
  zpos      = fgkSheight/2.0 - fgkSMpltT - 0.4 - kCOLwid/2.0; 
  parCOL[0] = kCOLwid   /2.0;
  parCOL[1] = fgkFlength/2.0;
  parCOL[2] = kCOLhgt   /2.0;
  TVirtualMC::GetMC()->Gsposp("UTC3",6+ 2*kNlayer,"UTF1", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsposp("UTC3",6+ 3*kNlayer,"UTF1",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsposp("UTC3",6+ 4*kNlayer,"UTF2", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsposp("UTC3",6+ 5*kNlayer,"UTF2",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  // In baby frame
  xpos      = fgkCwidth[5]/2.0 - kCOLhgt/2.0 - 3.1;
  ypos      = kBBSdz/2.0 - kBBMdz/2.0;
  zpos      = fgkSheight/2.0 - fgkSMpltT - 0.4 - kCOLwid/2.0; 
  parCOL[0] = kCOLwid/2.0;
  parCOL[1] = kBBSdz /2.0;
  parCOL[2] = kCOLhgt/2.0;
  TVirtualMC::GetMC()->Gsposp("UTC3",6+6*kNlayer,"BBTRD", xpos, ypos, zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsposp("UTC3",6+7*kNlayer,"BBTRD",-xpos, ypos, zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  // In back frame
  xpos      = fgkCwidth[5]/2.0 - kCOLhgt/2.0 - 1.3;
  ypos      = -kBFSdz/2.0 + kBFMdz/2.0;
  zpos      = fgkSheight/2.0 - fgkSMpltT - 0.4 - kCOLwid/2.0; 
  parCOL[0] = kCOLwid/2.0;
  parCOL[1] = kBFSdz /2.0;
  parCOL[2] = kCOLhgt/2.0;
  TVirtualMC::GetMC()->Gsposp("UTC3",6+6*kNlayer,"BFTRD", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);
  TVirtualMC::GetMC()->Gsposp("UTC3",6+7*kNlayer,"BFTRD",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parCOL,kNparCOL);

  //
  // The power bus bars
  //

  const Float_t kPWRwid  =  0.6;
  // Increase the height of the power bus bars to take into
  // account the material of additional cables, etc.
  const Float_t kPWRhgtA =  5.0 + 0.2;
  const Float_t kPWRhgtB =  5.0;
  const Float_t kPWRposx =  2.0;
  const Float_t kPWRposz =  0.1;
  const Int_t   kNparPWR =  3;
  Float_t parPWR[kNparPWR];
  parPWR[0] = 0.0;
  parPWR[1] = 0.0;
  parPWR[2] = 0.0;
  TVirtualMC::GetMC()->Gsvolu("UTP1","BOX ",idtmed[1325-1],parPWR,0);
  TVirtualMC::GetMC()->Gsvolu("UTP3","BOX ",idtmed[1325-1],parPWR,0);
  
  for (ilayer = 1; ilayer < kNlayer; ilayer++) { 

    // Along the chambers
    xpos      = fgkCwidth[ilayer]/2.0 + kPWRwid/2.0 + kPWRposx;
    ypos      = 0.0;
    zpos      = fgkVrocsm + fgkSMpltT - fgkCalZpos 
              + kPWRhgtA/2.0 - fgkSheight/2.0 + kPWRposz 
              + ilayer * (fgkCH + fgkVspace);
    parPWR[0] = kPWRwid   /2.0;
    parPWR[1] = fgkSlength/2.0;
    parPWR[2] = kPWRhgtA  /2.0;
    TVirtualMC::GetMC()->Gsposp("UTP1",ilayer           ,"UTI1", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parPWR,kNparPWR);
    TVirtualMC::GetMC()->Gsposp("UTP1",ilayer+   kNlayer,"UTI1",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parPWR,kNparPWR);
    TVirtualMC::GetMC()->Gsposp("UTP1",ilayer+ 6*kNlayer,"UTI2", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parPWR,kNparPWR);
    TVirtualMC::GetMC()->Gsposp("UTP1",ilayer+ 7*kNlayer,"UTI2",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parPWR,kNparPWR);
    TVirtualMC::GetMC()->Gsposp("UTP1",ilayer+ 8*kNlayer,"UTI3", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parPWR,kNparPWR);
    TVirtualMC::GetMC()->Gsposp("UTP1",ilayer+ 9*kNlayer,"UTI3",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parPWR,kNparPWR);
    TVirtualMC::GetMC()->Gsposp("UTP1",ilayer+10*kNlayer,"UTI4", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parPWR,kNparPWR);
    TVirtualMC::GetMC()->Gsposp("UTP1",ilayer+11*kNlayer,"UTI4",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parPWR,kNparPWR);

    // Front of supermodule
    xpos      = fgkCwidth[ilayer]/2.0 + kPWRwid/2.0 + kPWRposx;
    ypos      = 0.0;
    zpos      = fgkVrocsm + fgkSMpltT - fgkCalZpos 
              + kPWRhgtA/2.0 - fgkSheight/2.0 + kPWRposz 
              + ilayer * (fgkCH + fgkVspace);
    parPWR[0] = kPWRwid   /2.0;
    parPWR[1] = fgkFlength/2.0;
    parPWR[2] = kPWRhgtA  /2.0;
    TVirtualMC::GetMC()->Gsposp("UTP3",ilayer+ 2*kNlayer,"UTF1", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parPWR,kNparPWR);
    TVirtualMC::GetMC()->Gsposp("UTP3",ilayer+ 3*kNlayer,"UTF1",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parPWR,kNparPWR);
    TVirtualMC::GetMC()->Gsposp("UTP3",ilayer+ 4*kNlayer,"UTF2", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parPWR,kNparPWR);
    TVirtualMC::GetMC()->Gsposp("UTP3",ilayer+ 5*kNlayer,"UTF2",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parPWR,kNparPWR);

  }

  for (ilayer = 1; ilayer < kNlayer; ilayer++) { 

    // In baby frame
    xpos      = fgkCwidth[ilayer]/2.0 + kPWRwid/2.0 + kPWRposx - 2.5;
    ypos      = kBBSdz/2.0 - kBBMdz/2.0;
    zpos      = fgkVrocsm + fgkSMpltT - fgkCalZpos 
              + kPWRhgtB/2.0 - fgkSheight/2.0 + kPWRposz 
              + ilayer * (fgkCH + fgkVspace);
    parPWR[0] = kPWRwid /2.0;
    parPWR[1] = kBBSdz  /2.0;
    parPWR[2] = kPWRhgtB/2.0;
    TVirtualMC::GetMC()->Gsposp("UTP3",ilayer+6*kNlayer,"BBTRD", xpos, ypos, zpos
                      ,matrix[0],"ONLY",parPWR,kNparPWR);
    TVirtualMC::GetMC()->Gsposp("UTP3",ilayer+7*kNlayer,"BBTRD",-xpos, ypos, zpos
                      ,matrix[1],"ONLY",parPWR,kNparPWR);

  }

  for (ilayer = 1; ilayer < kNlayer; ilayer++) { 

    // In back frame
    xpos      = fgkCwidth[ilayer]/2.0 + kPWRwid/2.0 + kPWRposx - 0.3;
    ypos      = -kBFSdz/2.0 + kBFMdz/2.0;
    zpos      = fgkVrocsm + fgkSMpltT - fgkCalZpos 
              + kPWRhgtB/2.0 - fgkSheight/2.0 + kPWRposz 
              + ilayer * (fgkCH + fgkVspace);
    parPWR[0] = kPWRwid /2.0;
    parPWR[1] = kBFSdz  /2.0;
    parPWR[2] = kPWRhgtB/2.0;
    TVirtualMC::GetMC()->Gsposp("UTP3",ilayer+8*kNlayer,"BFTRD", xpos,ypos,zpos
                      ,matrix[0],"ONLY",parPWR,kNparPWR);
    TVirtualMC::GetMC()->Gsposp("UTP3",ilayer+9*kNlayer,"BFTRD",-xpos,ypos,zpos
                      ,matrix[1],"ONLY",parPWR,kNparPWR);

  }

  // The upper most layer
  // Along the chambers
  xpos      = fgkCwidth[5]/2.0 + kPWRhgtB/2.0 - 1.3;
  ypos      = 0.0;
  zpos      = fgkSheight/2.0 - fgkSMpltT - 0.6 - kPWRwid/2.0; 
  parPWR[0] = kPWRwid   /2.0;
  parPWR[1] = fgkSlength/2.0;
  parPWR[2] = kPWRhgtB  /2.0 ;
  TVirtualMC::GetMC()->Gsposp("UTP1",6           ,"UTI1", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  TVirtualMC::GetMC()->Gsposp("UTP1",6+   kNlayer,"UTI1",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  TVirtualMC::GetMC()->Gsposp("UTP1",6+ 6*kNlayer,"UTI2", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  TVirtualMC::GetMC()->Gsposp("UTP1",6+ 7*kNlayer,"UTI2",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  TVirtualMC::GetMC()->Gsposp("UTP1",6+ 8*kNlayer,"UTI3", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  TVirtualMC::GetMC()->Gsposp("UTP1",6+ 9*kNlayer,"UTI3",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  TVirtualMC::GetMC()->Gsposp("UTP1",6+10*kNlayer,"UTI4", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  TVirtualMC::GetMC()->Gsposp("UTP1",6+11*kNlayer,"UTI4",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  // Front of supermodules
  xpos      = fgkCwidth[5]/2.0 + kPWRhgtB/2.0 - 1.3;
  ypos      = 0.0;
  zpos      = fgkSheight/2.0 - fgkSMpltT - 0.6 - kPWRwid/2.0; 
  parPWR[0] = kPWRwid   /2.0;
  parPWR[1] = fgkFlength/2.0;
  parPWR[2] = kPWRhgtB  /2.0;
  TVirtualMC::GetMC()->Gsposp("UTP3",6+2*kNlayer,"UTF1", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  TVirtualMC::GetMC()->Gsposp("UTP3",6+3*kNlayer,"UTF1",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  TVirtualMC::GetMC()->Gsposp("UTP3",6+4*kNlayer,"UTF2", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  TVirtualMC::GetMC()->Gsposp("UTP3",6+5*kNlayer,"UTF2",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  // In baby frame
  xpos      = fgkCwidth[5]/2.0 + kPWRhgtB/2.0 - 3.0;
  ypos      = kBBSdz/2.0 - kBBMdz/2.0;
  zpos      = fgkSheight/2.0 - fgkSMpltT - 0.6 - kPWRwid/2.0; 
  parPWR[0] = kPWRwid /2.0;
  parPWR[1] = kBBSdz  /2.0;
  parPWR[2] = kPWRhgtB/2.0;
  TVirtualMC::GetMC()->Gsposp("UTP3",6+6*kNlayer,"BBTRD", xpos, ypos, zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  TVirtualMC::GetMC()->Gsposp("UTP3",6+7*kNlayer,"BBTRD",-xpos, ypos, zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  // In back frame
  xpos      = fgkCwidth[5]/2.0 + kPWRhgtB/2.0 - 1.3;
  ypos      = -kBFSdz/2.0 + kBFMdz/2.0;
  zpos      = fgkSheight/2.0 - fgkSMpltT - 0.6 - kPWRwid/2.0; 
  parPWR[0] = kPWRwid /2.0;
  parPWR[1] = kBFSdz  /2.0;
  parPWR[2] = kPWRhgtB/2.0;
  TVirtualMC::GetMC()->Gsposp("UTP3",6+8*kNlayer,"BFTRD", xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);
  TVirtualMC::GetMC()->Gsposp("UTP3",6+9*kNlayer,"BFTRD",-xpos,ypos,zpos
                    ,matrix[3],"ONLY",parPWR,kNparPWR);

  //
  // The gas tubes connecting the chambers in the super modules with holes
  // Material: Stainless steel
  //

  // PHOS holes
  parTube[0] = 0.0;
  parTube[1] = 2.2/2.0;
  parTube[2] = fgkClength[5][2]/2.0 - fgkHspace/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTG1","TUBE",idtmed[1308-1],parTube,kNparTube);
  parTube[0] = 0.0;
  parTube[1] = 2.1/2.0;
  parTube[2] = fgkClength[5][2]/2.0 - fgkHspace/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTG2","TUBE",idtmed[1309-1],parTube,kNparTube);
  xpos  = 0.0;
  ypos  = 0.0;
  zpos  = 0.0;
  TVirtualMC::GetMC()->Gspos("UTG2",1,"UTG1",xpos,ypos,zpos,0,"ONLY");
  for (ilayer = 0; ilayer < kNlayer; ilayer++) { 
    xpos      = fgkCwidth[ilayer]/2.0 + kCOLwid/2.0 - 1.5;
    ypos      = 0.0;
    zpos      = fgkVrocsm + fgkSMpltT + kCOLhgt/2.0 - fgkSheight/2.0 + 5.0 
              + ilayer * (fgkCH + fgkVspace);
    TVirtualMC::GetMC()->Gspos("UTG1",1+ilayer,"UTI3", xpos, ypos, zpos,matrix[4],"ONLY");
    TVirtualMC::GetMC()->Gspos("UTG1",7+ilayer,"UTI3",-xpos, ypos, zpos,matrix[4],"ONLY");
  }
  // Missing L4S4 chamber in sector 17
  parTube[0] = 0.0;
  parTube[1] = 2.2/2.0;
  parTube[2] = fgkClength[4][4]/2.0 - fgkHspace/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTG3","TUBE",idtmed[1308-1],parTube,kNparTube);
  parTube[0] = 0.0;
  parTube[1] = 2.1/2.0;
  parTube[2] = fgkClength[4][4]/2.0 - fgkHspace/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTG4","TUBE",idtmed[1309-1],parTube,kNparTube);
  xpos  = 0.0;
  ypos  = 0.0;
  zpos  = 0.0;
  TVirtualMC::GetMC()->Gspos("UTG4",1,"UTG3",xpos,ypos,zpos,0,"ONLY");
  xpos  = fgkCwidth[4]/2.0 + kCOLwid/2.0 - 1.5;
  ypos  = -fgkClength[4][0]/2.0 - fgkClength[4][1] - fgkClength[4][2]/2.0;
  zpos  = fgkVrocsm + fgkSMpltT + kCOLhgt/2.0 - fgkSheight/2.0 + 5.0 + 4 * (fgkCH + fgkVspace);
  TVirtualMC::GetMC()->Gspos("UTG3",1,"UTI4", xpos, ypos, zpos,matrix[4],"ONLY");
  TVirtualMC::GetMC()->Gspos("UTG4",2,"UTI4",-xpos, ypos, zpos,matrix[4],"ONLY");

  //
  // The volumes for the services at the chambers
  //

  const Int_t kNparServ = 3;
  Float_t parServ[kNparServ];

  for (istack = 0; istack < kNstack; istack++) {
    for (ilayer = 0; ilayer < kNlayer; ilayer++) {

      Int_t iDet = GetDetectorSec(ilayer,istack);

      snprintf(cTagV,kTag,"UU%02d",iDet);
      parServ[0] = fgkCwidth[ilayer]         /2.0;
      parServ[1] = fgkClength[ilayer][istack]/2.0 - fgkHspace/2.0;
      parServ[2] = fgkCsvH                 /2.0;
      TVirtualMC::GetMC()->Gsvolu(cTagV,"BOX",idtmed[1302-1],parServ,kNparServ);

    }
  }

  //
  // The cooling pipes inside the service volumes
  //

  // The cooling pipes
  parTube[0] =  0.0;
  parTube[1] =  0.0;
  parTube[2] =  0.0;
  TVirtualMC::GetMC()->Gsvolu("UTCP","TUBE",idtmed[1324-1],parTube,0);
  // The cooling water
  parTube[0] =  0.0;
  parTube[1] =  0.2/2.0;
  parTube[2] = -1.0;
  TVirtualMC::GetMC()->Gsvolu("UTCH","TUBE",idtmed[1314-1],parTube,kNparTube);
  // Water inside the cooling pipe
  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  TVirtualMC::GetMC()->Gspos("UTCH",1,"UTCP",xpos,ypos,zpos,0,"ONLY");

  // Position the cooling pipes in the mother volume
  for (istack = 0; istack < kNstack; istack++) {
    for (ilayer = 0; ilayer < kNlayer; ilayer++) {
      Int_t   iDet    = GetDetectorSec(ilayer,istack);
      Int_t   iCopy   = GetDetector(ilayer,istack,0) * 100;
      Int_t   nMCMrow = GetRowMax(ilayer,istack,0);
      Float_t ySize   = (GetChamberLength(ilayer,istack) - 2.0*fgkRpadW) 
                      / ((Float_t) nMCMrow);
      snprintf(cTagV,kTag,"UU%02d",iDet);
      for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
        xpos   = 0.0;
        ypos   = (0.5 + iMCMrow) * ySize 
               - fgkClength[ilayer][istack]/2.0 + fgkHspace/2.0;
        zpos   = 0.0 + 0.742/2.0;                 
	// The cooling pipes
        parTube[0] = 0.0;
        parTube[1] = 0.3/2.0; // Thickness of the cooling pipes
        parTube[2] = fgkCwidth[ilayer]/2.0;
        TVirtualMC::GetMC()->Gsposp("UTCP",iCopy+iMCMrow,cTagV,xpos,ypos,zpos
                          ,matrix[2],"ONLY",parTube,kNparTube);
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
  TVirtualMC::GetMC()->Gsvolu("UTPL","TUBE",idtmed[1305-1],parTube,0);

  // Position the power lines in the mother volume
  for (istack = 0; istack < kNstack; istack++) {
    for (ilayer = 0; ilayer < kNlayer; ilayer++) {
      Int_t   iDet    = GetDetectorSec(ilayer,istack);
      Int_t   iCopy   = GetDetector(ilayer,istack,0) * 100;
      Int_t   nMCMrow = GetRowMax(ilayer,istack,0);
      Float_t ySize   = (GetChamberLength(ilayer,istack) - 2.0*fgkRpadW) 
                      / ((Float_t) nMCMrow);
      snprintf(cTagV,kTag,"UU%02d",iDet);
      for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
        xpos       = 0.0;
        ypos       = (0.5 + iMCMrow) * ySize - 1.0 
                   - fgkClength[ilayer][istack]/2.0 + fgkHspace/2.0;
        zpos       = -0.4 + 0.742/2.0;
        parTube[0] = 0.0;
        parTube[1] = 0.2/2.0; // Thickness of the power lines
        parTube[2] = fgkCwidth[ilayer]/2.0;
        TVirtualMC::GetMC()->Gsposp("UTPL",iCopy+iMCMrow,cTagV,xpos,ypos,zpos
                          ,matrix[2],"ONLY",parTube,kNparTube);
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
  const Float_t kMCMcuTh = 0.0025;
  const Float_t kMCMsiTh = 0.03;
  const Float_t kMCMcoTh = 0.04;

  // The mother volume for the MCMs (air)
  const Int_t kNparMCM = 3;
  Float_t parMCM[kNparMCM];
  parMCM[0] = kMCMx   /2.0;
  parMCM[1] = kMCMy   /2.0;
  parMCM[2] = kMCMz   /2.0;
  TVirtualMC::GetMC()->Gsvolu("UMCM","BOX",idtmed[1302-1],parMCM,kNparMCM);

  // The MCM carrier G10 layer
  parMCM[0] = kMCMx   /2.0;
  parMCM[1] = kMCMy   /2.0;
  parMCM[2] = kMCMpcTh/2.0;
  TVirtualMC::GetMC()->Gsvolu("UMC1","BOX",idtmed[1319-1],parMCM,kNparMCM);
  // The MCM carrier Cu layer
  parMCM[0] = kMCMx   /2.0;
  parMCM[1] = kMCMy   /2.0;
  parMCM[2] = kMCMcuTh/2.0;
  TVirtualMC::GetMC()->Gsvolu("UMC2","BOX",idtmed[1318-1],parMCM,kNparMCM);
  // The silicon of the chips
  parMCM[0] = kMCMx   /2.0;
  parMCM[1] = kMCMy   /2.0;
  parMCM[2] = kMCMsiTh/2.0;
  TVirtualMC::GetMC()->Gsvolu("UMC3","BOX",idtmed[1320-1],parMCM,kNparMCM);
  // The aluminum of the cooling plates
  parMCM[0] = kMCMx   /2.0;
  parMCM[1] = kMCMy   /2.0;
  parMCM[2] = kMCMcoTh/2.0;
  TVirtualMC::GetMC()->Gsvolu("UMC4","BOX",idtmed[1324-1],parMCM,kNparMCM);

  // Put the MCM material inside the MCM mother volume
  xpos  =  0.0;
  ypos  =  0.0;
  zpos  = -kMCMz   /2.0 + kMCMpcTh/2.0;
  TVirtualMC::GetMC()->Gspos("UMC1",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  kMCMpcTh/2.0 + kMCMcuTh/2.0;
  TVirtualMC::GetMC()->Gspos("UMC2",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  kMCMcuTh/2.0 + kMCMsiTh/2.0;
  TVirtualMC::GetMC()->Gspos("UMC3",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  kMCMsiTh/2.0 + kMCMcoTh/2.0;
  TVirtualMC::GetMC()->Gspos("UMC4",1,"UMCM",xpos,ypos,zpos,0,"ONLY");

  // Position the MCMs in the mother volume
  for (istack = 0; istack < kNstack; istack++) {
    for (ilayer = 0; ilayer < kNlayer; ilayer++) {
      Int_t   iDet    = GetDetectorSec(ilayer,istack);
      Int_t   iCopy   = GetDetector(ilayer,istack,0) * 1000;
      Int_t   nMCMrow = GetRowMax(ilayer,istack,0);
      Float_t ySize   = (GetChamberLength(ilayer,istack) - 2.0*fgkRpadW)
                      / ((Float_t) nMCMrow);
      Int_t   nMCMcol = 8;
      Float_t xSize   = (GetChamberWidth(ilayer)         - 2.0*fgkCpadW)
                      / ((Float_t) nMCMcol + 6);             // Introduce 6 gaps
      Int_t   iMCM[8] = {  1,  2,  3,  5,  8,  9, 10, 12 };  // 0..7 MCM + 6 gap structure
      snprintf(cTagV,kTag,"UU%02d",iDet);
      for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
        for (Int_t iMCMcol = 0; iMCMcol < nMCMcol; iMCMcol++) {
          xpos      = (0.5 + iMCM[iMCMcol]) * xSize + 1.0
                    - fgkCwidth[ilayer]/2.0;
          ypos      = (0.5 + iMCMrow) * ySize + 1.0
                    - fgkClength[ilayer][istack]/2.0 + fgkHspace/2.0;
          zpos      = -0.4 + 0.742/2.0;
          TVirtualMC::GetMC()->Gspos("UMCM",iCopy+iMCMrow*10+iMCMcol,cTagV
                           ,xpos,ypos,zpos,0,"ONLY");
	  // Add two additional smaller cooling pipes on top of the MCMs
	  // to mimic the meandering structure
          xpos      = (0.5 + iMCM[iMCMcol]) * xSize + 1.0
                    - fgkCwidth[ilayer]/2.0;
          ypos      = (0.5 + iMCMrow) * ySize
                    - fgkClength[ilayer][istack]/2.0 + fgkHspace/2.0;
          zpos      = 0.0 + 0.742/2.0;                 
          parTube[0] = 0.0;
          parTube[1] = 0.3/2.0; // Thickness of the cooling pipes
          parTube[2] = kMCMx/2.0;
          TVirtualMC::GetMC()->Gsposp("UTCP",iCopy+iMCMrow*10+iMCMcol+ 50,cTagV
                            ,xpos,ypos+1.0,zpos
                            ,matrix[2],"ONLY",parTube,kNparTube);
          TVirtualMC::GetMC()->Gsposp("UTCP",iCopy+iMCMrow*10+iMCMcol+500,cTagV
                            ,xpos,ypos+2.0,zpos
                            ,matrix[2],"ONLY",parTube,kNparTube);

        }
      }

    }
  }

  //
  // The DCS boards
  //

  const Float_t kDCSx    =  9.0;
  const Float_t kDCSy    = 14.5;
  const Float_t kDCSz    =  0.3;
  
  const Float_t kDCSpcTh =  0.15;
  const Float_t kDCScuTh =  0.01;
  const Float_t kDCScoTh =  0.04;

  // The mother volume for the DCSs (air)
  const Int_t kNparDCS = 3;
  Float_t parDCS[kNparDCS];
  parDCS[0] = kDCSx   /2.0;
  parDCS[1] = kDCSy   /2.0;
  parDCS[2] = kDCSz   /2.0;
  TVirtualMC::GetMC()->Gsvolu("UDCS","BOX",idtmed[1302-1],parDCS,kNparDCS);

  // The DCS carrier G10 layer
  parDCS[0] = kDCSx   /2.0;
  parDCS[1] = kDCSy   /2.0;
  parDCS[2] = kDCSpcTh/2.0;
  TVirtualMC::GetMC()->Gsvolu("UDC1","BOX",idtmed[1319-1],parDCS,kNparDCS);
  // The DCS carrier Cu layer
  parDCS[0] = kDCSx   /2.0;
  parDCS[1] = kDCSy   /2.0;
  parDCS[2] = kDCScuTh/2.0;
  TVirtualMC::GetMC()->Gsvolu("UDC2","BOX",idtmed[1318-1],parDCS,kNparDCS);
  // The aluminum of the cooling plates
  parDCS[0] = 5.0     /2.0;
  parDCS[1] = 5.0     /2.0;
  parDCS[2] = kDCScoTh/2.0;
  TVirtualMC::GetMC()->Gsvolu("UDC3","BOX",idtmed[1324-1],parDCS,kNparDCS);

  // Put the DCS material inside the DCS mother volume
  xpos  =  0.0;
  ypos  =  0.0;
  zpos  = -kDCSz   /2.0 + kDCSpcTh/2.0;
  TVirtualMC::GetMC()->Gspos("UDC1",1,"UDCS",xpos,ypos,zpos,0,"ONLY");
  zpos +=  kDCSpcTh/2.0 + kDCScuTh/2.0;
  TVirtualMC::GetMC()->Gspos("UDC2",1,"UDCS",xpos,ypos,zpos,0,"ONLY");
  zpos +=  kDCScuTh/2.0 + kDCScoTh/2.0;
  TVirtualMC::GetMC()->Gspos("UDC3",1,"UDCS",xpos,ypos,zpos,0,"ONLY");

  // Put the DCS board in the chamber services mother volume
  for (istack = 0; istack < kNstack; istack++) {
    for (ilayer = 0; ilayer < kNlayer; ilayer++) {
      Int_t   iDet    = GetDetectorSec(ilayer,istack);
      Int_t   iCopy   = iDet + 1;
      xpos =  fgkCwidth[ilayer]/2.0 - 1.9 * (GetChamberLength(ilayer,istack) - 2.0*fgkRpadW) 
                                        / ((Float_t) GetRowMax(ilayer,istack,0));
      ypos =  0.05 * fgkClength[ilayer][istack];
      zpos =  kDCSz/2.0 - fgkCsvH/2.0;
      snprintf(cTagV,kTag,"UU%02d",iDet);
      TVirtualMC::GetMC()->Gspos("UDCS",iCopy,cTagV,xpos,ypos,zpos,0,"ONLY");
    }
  }

  //
  // The ORI boards
  //

  const Float_t kORIx    =  4.2;
  const Float_t kORIy    = 13.5;
  const Float_t kORIz    =  0.3;
  
  const Float_t kORIpcTh =  0.15;
  const Float_t kORIcuTh =  0.01;
  const Float_t kORIcoTh =  0.04;

  // The mother volume for the ORIs (air)
  const Int_t kNparORI = 3;
  Float_t parORI[kNparORI];
  parORI[0] = kORIx   /2.0;
  parORI[1] = kORIy   /2.0;
  parORI[2] = kORIz   /2.0;
  TVirtualMC::GetMC()->Gsvolu("UORI","BOX",idtmed[1302-1],parORI,kNparORI);

  // The ORI carrier G10 layer
  parORI[0] = kORIx   /2.0;
  parORI[1] = kORIy   /2.0;
  parORI[2] = kORIpcTh/2.0;
  TVirtualMC::GetMC()->Gsvolu("UOR1","BOX",idtmed[1319-1],parORI,kNparORI);
  // The ORI carrier Cu layer
  parORI[0] = kORIx   /2.0;
  parORI[1] = kORIy   /2.0;
  parORI[2] = kORIcuTh/2.0;
  TVirtualMC::GetMC()->Gsvolu("UOR2","BOX",idtmed[1318-1],parORI,kNparORI);
  // The aluminum of the cooling plates
  parORI[0] = kORIx   /2.0;
  parORI[1] = kORIy   /2.0;
  parORI[2] = kORIcoTh/2.0;
  TVirtualMC::GetMC()->Gsvolu("UOR3","BOX",idtmed[1324-1],parORI,kNparORI);

  // Put the ORI material inside the ORI mother volume
  xpos  =  0.0;
  ypos  =  0.0;
  zpos  = -kORIz   /2.0 + kORIpcTh/2.0;
  TVirtualMC::GetMC()->Gspos("UOR1",1,"UORI",xpos,ypos,zpos,0,"ONLY");
  zpos +=  kORIpcTh/2.0 + kORIcuTh/2.0;
  TVirtualMC::GetMC()->Gspos("UOR2",1,"UORI",xpos,ypos,zpos,0,"ONLY");
  zpos +=  kORIcuTh/2.0 + kORIcoTh/2.0;
  TVirtualMC::GetMC()->Gspos("UOR3",1,"UORI",xpos,ypos,zpos,0,"ONLY");

  // Put the ORI board in the chamber services mother volume
  for (istack = 0; istack < kNstack; istack++) {
    for (ilayer = 0; ilayer < kNlayer; ilayer++) {
      Int_t   iDet    = GetDetectorSec(ilayer,istack);
      Int_t   iCopy   = iDet + 1;
      xpos =  fgkCwidth[ilayer]/2.0 - 1.92 * (GetChamberLength(ilayer,istack) - 2.0*fgkRpadW) 
                                        / ((Float_t) GetRowMax(ilayer,istack,0));
      ypos = -16.0;
      zpos =  kORIz/2.0 - fgkCsvH/2.0;
      snprintf(cTagV,kTag,"UU%02d",iDet);
      TVirtualMC::GetMC()->Gspos("UORI",iCopy      ,cTagV,xpos,ypos,zpos,0,"ONLY");
      xpos = -fgkCwidth[ilayer]/2.0 + 3.8 * (GetChamberLength(ilayer,istack) - 2.0*fgkRpadW) 
                                        / ((Float_t) GetRowMax(ilayer,istack,0));
      ypos = -16.0;
      zpos =  kORIz/2.0 - fgkCsvH/2.0;
      snprintf(cTagV,kTag,"UU%02d",iDet);
      TVirtualMC::GetMC()->Gspos("UORI",iCopy+kNdet,cTagV,xpos,ypos,zpos,0,"ONLY");
    }
  }

  //
  // Services in front of the super module
  //

  // Gas in-/outlet pipes (INOX)
  parTube[0] = 0.0;
  parTube[1] = 0.0;
  parTube[2] = 0.0;
  TVirtualMC::GetMC()->Gsvolu("UTG3","TUBE",idtmed[1308-1],parTube,0);
  // The gas inside the in-/outlet pipes (Xe)
  parTube[0] =  0.0;
  parTube[1] =  1.2/2.0;
  parTube[2] = -1.0;
  TVirtualMC::GetMC()->Gsvolu("UTG4","TUBE",idtmed[1309-1],parTube,kNparTube);
  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  TVirtualMC::GetMC()->Gspos("UTG4",1,"UTG3",xpos,ypos,zpos,0,"ONLY");
  for (ilayer = 0; ilayer < kNlayer-1; ilayer++) { 
    xpos       = 0.0;
    ypos       = fgkClength[ilayer][2]/2.0 
               + fgkClength[ilayer][1] 
               + fgkClength[ilayer][0];
    zpos       = 9.0 - fgkSheight/2.0
               + ilayer * (fgkCH + fgkVspace);
    parTube[0] = 0.0;
    parTube[1] = 1.5/2.0;
    parTube[2] = fgkCwidth[ilayer]/2.0 - 2.5;
    TVirtualMC::GetMC()->Gsposp("UTG3",ilayer+1          ,"UTI1", xpos, ypos, zpos
                      ,matrix[2],"ONLY",parTube,kNparTube);
    TVirtualMC::GetMC()->Gsposp("UTG3",ilayer+1+1*kNlayer,"UTI1", xpos,-ypos, zpos
                      ,matrix[2],"ONLY",parTube,kNparTube);
    TVirtualMC::GetMC()->Gsposp("UTG3",ilayer+1+2*kNlayer,"UTI2", xpos, ypos, zpos
                      ,matrix[2],"ONLY",parTube,kNparTube);
    TVirtualMC::GetMC()->Gsposp("UTG3",ilayer+1+3*kNlayer,"UTI2", xpos,-ypos, zpos
                      ,matrix[2],"ONLY",parTube,kNparTube);
    TVirtualMC::GetMC()->Gsposp("UTG3",ilayer+1+4*kNlayer,"UTI3", xpos, ypos, zpos
                      ,matrix[2],"ONLY",parTube,kNparTube);
    TVirtualMC::GetMC()->Gsposp("UTG3",ilayer+1+5*kNlayer,"UTI3", xpos,-ypos, zpos
                      ,matrix[2],"ONLY",parTube,kNparTube);
    TVirtualMC::GetMC()->Gsposp("UTG3",ilayer+1+6*kNlayer,"UTI4", xpos, ypos, zpos
                      ,matrix[2],"ONLY",parTube,kNparTube);
    TVirtualMC::GetMC()->Gsposp("UTG3",ilayer+1+7*kNlayer,"UTI4", xpos,-ypos, zpos
                      ,matrix[2],"ONLY",parTube,kNparTube);
  }

  // Gas distribution box
  parBox[0] = 14.50/2.0;
  parBox[1] =  4.52/2.0;
  parBox[2] =  5.00/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTGD","BOX ",idtmed[1308-1],parBox,kNparBox);
  parBox[0] = 14.50/2.0;
  parBox[1] =  4.00/2.0;
  parBox[2] =  4.40/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTGI","BOX ",idtmed[1309-1],parBox,kNparBox);
  parTube[0] = 0.0;
  parTube[1] = 4.0/2.0;
  parTube[2] = 8.0/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTGT","TUBE",idtmed[1308-1],parTube,kNparTube);
  parTube[0] = 0.0;
  parTube[1] = 3.4/2.0;
  parTube[2] = 8.0/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTGG","TUBE",idtmed[1309-1],parTube,kNparTube);
  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  TVirtualMC::GetMC()->Gspos("UTGI",1,"UTGD",xpos,ypos,zpos,        0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTGG",1,"UTGT",xpos,ypos,zpos,        0,"ONLY");
  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  TVirtualMC::GetMC()->Gspos("UTGD",1,"UTF1",xpos,ypos,zpos,        0,"ONLY");
  xpos =  -3.0;
  ypos =   0.0;
  zpos =   6.5;
  TVirtualMC::GetMC()->Gspos("UTGT",1,"UTF1",xpos,ypos,zpos,        0,"ONLY");
  xpos = -11.25;
  ypos =   0.0;
  zpos =   0.5;
  TVirtualMC::GetMC()->Gspos("UTGT",3,"UTF1",xpos,ypos,zpos,matrix[2],"ONLY");
  xpos =  11.25;
  ypos =   0.0;
  zpos =   0.5;
  TVirtualMC::GetMC()->Gspos("UTGT",5,"UTF1",xpos,ypos,zpos,matrix[2],"ONLY");

  // Cooling manifolds
  parBox[0]  =  5.0/2.0;
  parBox[1]  = 23.0/2.0;
  parBox[2]  = 70.0/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTCM","BOX ",idtmed[1302-1],parBox,kNparBox);
  parBox[0]  =  5.0/2.0;
  parBox[1]  =  5.0/2.0;
  parBox[2]  = 70.0/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTCA","BOX ",idtmed[1308-1],parBox,kNparBox);
  parBox[0]  =  5.0/2.0 - 0.3;
  parBox[1]  =  5.0/2.0 - 0.3;
  parBox[2]  = 70.0/2.0 - 0.3;
  TVirtualMC::GetMC()->Gsvolu("UTCW","BOX ",idtmed[1314-1],parBox,kNparBox);
  xpos       =  0.0;
  ypos       =  0.0;
  zpos       =  0.0;
  TVirtualMC::GetMC()->Gspos("UTCW",1,"UTCA", xpos, ypos, zpos,        0,"ONLY");
  xpos       =  0.0;
  ypos       =  5.0/2.0 - 23.0/2.0;
  zpos       =  0.0;
  TVirtualMC::GetMC()->Gspos("UTCA",1,"UTCM", xpos, ypos, zpos,        0,"ONLY");
  parTube[0] =  0.0;
  parTube[1] =  3.0/2.0;
  parTube[2] = 18.0/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTCO","TUBE",idtmed[1308-1],parTube,kNparTube);
  parTube[0] =  0.0;
  parTube[1] =  3.0/2.0 - 0.3;
  parTube[2] = 18.0/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTCL","TUBE",idtmed[1314-1],parTube,kNparTube);
  xpos       =  0.0;
  ypos       =  0.0;
  zpos       =  0.0;
  TVirtualMC::GetMC()->Gspos("UTCL",1,"UTCO", xpos, ypos, zpos,        0,"ONLY");
  xpos       =  0.0;
  ypos       =  2.5;
  zpos       = -70.0/2.0 + 7.0;
  TVirtualMC::GetMC()->Gspos("UTCO",1,"UTCM", xpos, ypos, zpos,matrix[4],"ONLY");
  zpos      +=  7.0;
  TVirtualMC::GetMC()->Gspos("UTCO",2,"UTCM", xpos, ypos, zpos,matrix[4],"ONLY");
  zpos      +=  7.0;
  TVirtualMC::GetMC()->Gspos("UTCO",3,"UTCM", xpos, ypos, zpos,matrix[4],"ONLY");
  zpos      +=  7.0;
  TVirtualMC::GetMC()->Gspos("UTCO",4,"UTCM", xpos, ypos, zpos,matrix[4],"ONLY");
  zpos      +=  7.0;
  TVirtualMC::GetMC()->Gspos("UTCO",5,"UTCM", xpos, ypos, zpos,matrix[4],"ONLY");
  zpos      +=  7.0;
  TVirtualMC::GetMC()->Gspos("UTCO",6,"UTCM", xpos, ypos, zpos,matrix[4],"ONLY");
  zpos      +=  7.0;
  TVirtualMC::GetMC()->Gspos("UTCO",7,"UTCM", xpos, ypos, zpos,matrix[4],"ONLY");
  zpos      +=  7.0;
  TVirtualMC::GetMC()->Gspos("UTCO",8,"UTCM", xpos, ypos, zpos,matrix[4],"ONLY");

  xpos = 40.0;
  ypos =  fgkFlength/2.0 - 23.0/2.0;
  zpos =  0.0;
  TVirtualMC::GetMC()->Gspos("UTCM",1,"UTF1", xpos, ypos, zpos,matrix[0],"ONLY");
  TVirtualMC::GetMC()->Gspos("UTCM",2,"UTF1",-xpos, ypos, zpos,matrix[1],"ONLY");
  TVirtualMC::GetMC()->Gspos("UTCM",3,"UTF2", xpos,-ypos, zpos,matrix[5],"ONLY");
  TVirtualMC::GetMC()->Gspos("UTCM",4,"UTF2",-xpos,-ypos, zpos,matrix[6],"ONLY");

  // Power connection boards (Cu)
  parBox[0] =  0.5/2.0;
  parBox[1] = 15.0/2.0;
  parBox[2] =  7.0/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTPC","BOX ",idtmed[1325-1],parBox,kNparBox);
  for (ilayer = 0; ilayer < kNlayer-1; ilayer++) { 
    xpos      = fgkCwidth[ilayer]/2.0 + kPWRwid/2.0;
    ypos      = 0.0;
    zpos      = fgkVrocsm + fgkSMpltT + kPWRhgtA/2.0 - fgkSheight/2.0 + kPWRposz 
              + (ilayer+1) * (fgkCH + fgkVspace);
    TVirtualMC::GetMC()->Gspos("UTPC",ilayer        ,"UTF1", xpos,ypos,zpos,matrix[0],"ONLY");
    TVirtualMC::GetMC()->Gspos("UTPC",ilayer+kNlayer,"UTF1",-xpos,ypos,zpos,matrix[1],"ONLY");
  }
  xpos      = fgkCwidth[5]/2.0 + kPWRhgtA/2.0 - 2.0;
  ypos      = 0.0;
  zpos      = fgkSheight/2.0 - fgkSMpltT - 2.0; 
  TVirtualMC::GetMC()->Gspos("UTPC",5        ,"UTF1", xpos,ypos,zpos,matrix[3],"ONLY");
  TVirtualMC::GetMC()->Gspos("UTPC",5+kNlayer,"UTF1",-xpos,ypos,zpos,matrix[3],"ONLY");

  // Power connection panel (Al)
  parBox[0] = 60.0/2.0;
  parBox[1] = 10.0/2.0;
  parBox[2] =  3.0/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTPP","BOX ",idtmed[1301-1],parBox,kNparBox);
  xpos      =  0.0;
  ypos      =  0.0;
  zpos      = 18.0;
  TVirtualMC::GetMC()->Gspos("UTPP",1,"UTF1", xpos,ypos,zpos,0,"ONLY");

  //
  // Electronics boxes
  //

  // Casing (INOX)
  parBox[0] = 60.0/2.0;
  parBox[1] = 10.0/2.0;
  parBox[2] =  6.0/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTE1","BOX ",idtmed[1308-1],parBox,kNparBox);
  // Interior (air)
  parBox[0] = parBox[0] - 0.5;
  parBox[1] = parBox[1] - 0.5;
  parBox[2] = parBox[2] - 0.5;
  TVirtualMC::GetMC()->Gsvolu("UTE2","BOX ",idtmed[1302-1],parBox,kNparBox);
  xpos      = 0.0;
  ypos      = 0.0;
  zpos      = 0.0;
  TVirtualMC::GetMC()->Gspos("UTE2",1,"UTE1",xpos,ypos,zpos,0,"ONLY");
  xpos      = 0.0;
  ypos      =  fgkSlength/2.0 - 10.0/2.0 - 3.0;
  zpos      = -fgkSheight/2.0 +  6.0/2.0 + 1.0;
  TVirtualMC::GetMC()->Gspos("UTE1",1,"UTI1", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTE1",2,"UTI2", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTE1",3,"UTI3", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTE1",4,"UTI4", xpos,ypos,zpos,0,"ONLY");

  // Casing (INOX)
  parBox[0] = 50.0/2.0;
  parBox[1] = 15.0/2.0;
  parBox[2] = 20.0/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTE3","BOX ",idtmed[1308-1],parBox,kNparBox);
  // Interior (air)
  parBox[0] = parBox[0] - 0.5;
  parBox[1] = parBox[1] - 0.5;
  parBox[2] = parBox[2] - 0.5;
  TVirtualMC::GetMC()->Gsvolu("UTE4","BOX ",idtmed[1302-1],parBox,kNparBox);
  xpos      = 0.0;
  ypos      = 0.0;
  zpos      = 0.0;
  TVirtualMC::GetMC()->Gspos("UTE4",1,"UTE3",xpos,ypos,zpos,0,"ONLY");
  xpos      = 0.0;
  ypos      = -fgkSlength/2.0 + 15.0/2.0 + 3.0;
  zpos      = -fgkSheight/2.0 + 20.0/2.0 + 1.0;
  TVirtualMC::GetMC()->Gspos("UTE3",1,"UTI1", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTE3",2,"UTI2", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTE3",3,"UTI3", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTE3",4,"UTI4", xpos,ypos,zpos,0,"ONLY");

  // Casing (INOX)
  parBox[0] = 20.0/2.0;
  parBox[1] =  7.0/2.0;
  parBox[2] = 20.0/2.0;
  TVirtualMC::GetMC()->Gsvolu("UTE5","BOX ",idtmed[1308-1],parBox,kNparBox);
  // Interior (air)
  parBox[0] = parBox[0] - 0.5;
  parBox[1] = parBox[1] - 0.5;
  parBox[2] = parBox[2] - 0.5;
  TVirtualMC::GetMC()->Gsvolu("UTE6","BOX ",idtmed[1302-1],parBox,kNparBox);
  xpos      = 0.0;
  ypos      = 0.0;
  zpos      = 0.0;
  TVirtualMC::GetMC()->Gspos("UTE6",1,"UTE5",xpos,ypos,zpos,0,"ONLY");
  xpos      = 20.0;
  ypos      = -fgkSlength/2.0 +  7.0/2.0 + 3.0;
  zpos      = 0.0;
  TVirtualMC::GetMC()->Gspos("UTE5",1,"UTI1", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTE5",2,"UTI2", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTE5",3,"UTI3", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTE5",4,"UTI4", xpos,ypos,zpos,0,"ONLY");
  xpos      = -xpos;
  TVirtualMC::GetMC()->Gspos("UTE5",5,"UTI1", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTE5",6,"UTI2", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTE5",7,"UTI3", xpos,ypos,zpos,0,"ONLY");
  TVirtualMC::GetMC()->Gspos("UTE5",8,"UTI4", xpos,ypos,zpos,0,"ONLY");

}

//_____________________________________________________________________________
void AliTRDgeometry::AssembleChamber(Int_t ilayer, Int_t istack)
{
  //
  // Group volumes UA, UD, UF, UU into an assembly that defines the
  // alignable volume of a single readout chamber
  //

  const Int_t kTag = 100;
  Char_t  cTagV[kTag];
  Char_t  cTagM[kTag];

  Double_t xpos = 0.0;
  Double_t ypos = 0.0;
  Double_t zpos = 0.0;

  Int_t idet = GetDetectorSec(ilayer,istack);

  // Create the assembly for a given ROC
  snprintf(cTagM,kTag,"UT%02d",idet);
  TGeoVolume *roc = new TGeoVolumeAssembly(cTagM);

  // Add the lower part of the chamber (aluminum frame),
  // including radiator and drift region
  xpos = 0.0;
  ypos = 0.0;
  zpos = fgkCraH/2.0 + fgkCdrH/2.0 - fgkCHsv/2.0;
  snprintf(cTagV,kTag,"UA%02d",idet);
  TGeoVolume *rocA = gGeoManager->GetVolume(cTagV);
  roc->AddNode(rocA,1,new TGeoTranslation(xpos,ypos,zpos));

  // Add the additional aluminum ledges
  xpos = fgkCwidth[ilayer]/2.0 + fgkCalWmod/2.0;
  ypos = 0.0;
  zpos = fgkCraH + fgkCdrH - fgkCalZpos - fgkCalHmod/2.0 - fgkCHsv/2.0;
  snprintf(cTagV,kTag,"UZ%02d",idet);
  TGeoVolume *rocZ = gGeoManager->GetVolume(cTagV);
  roc->AddNode(rocZ,1,new TGeoTranslation( xpos,ypos,zpos));
  roc->AddNode(rocZ,2,new TGeoTranslation(-xpos,ypos,zpos));

  // Add the additional wacosit ledges
  xpos = fgkCwidth[ilayer]/2.0 + fgkCwsW/2.0;
  ypos = 0.0;
  zpos = fgkCraH + fgkCdrH - fgkCwsH/2.0 - fgkCHsv/2.0;
  snprintf(cTagV,kTag,"UP%02d",idet);
  TGeoVolume *rocP = gGeoManager->GetVolume(cTagV);
  roc->AddNode(rocP,1,new TGeoTranslation( xpos,ypos,zpos));
  roc->AddNode(rocP,2,new TGeoTranslation(-xpos,ypos,zpos));

  // Add the middle part of the chamber (G10 frame),
  // including amplification region
  xpos = 0.0;
  ypos = 0.0;
  zpos = fgkCamH/2.0 + fgkCraH + fgkCdrH - fgkCHsv/2.0;
  snprintf(cTagV,kTag,"UD%02d",idet);
  TGeoVolume *rocD = gGeoManager->GetVolume(cTagV);
  roc->AddNode(rocD,1,new TGeoTranslation(xpos,ypos,zpos));

  // Add the upper part of the chamber (aluminum frame),
  // including back panel and FEE
  xpos = 0.0;
  ypos = 0.0;
  zpos = fgkCroH/2.0 + fgkCamH + fgkCraH + fgkCdrH - fgkCHsv/2.0;
  snprintf(cTagV,kTag,"UF%02d",idet);
  TGeoVolume *rocF = gGeoManager->GetVolume(cTagV);
  roc->AddNode(rocF,1,new TGeoTranslation(xpos,ypos,zpos));

  // Add the volume with services on top of the back panel
  xpos = 0.0;
  ypos = 0.0;
  zpos = fgkCsvH/2.0 + fgkCroH + fgkCamH + fgkCraH + fgkCdrH - fgkCHsv/2.0;
  snprintf(cTagV,kTag,"UU%02d",idet);
  TGeoVolume *rocU = gGeoManager->GetVolume(cTagV);
  roc->AddNode(rocU,1,new TGeoTranslation(xpos,ypos,zpos));

  // Place the ROC assembly into the super modules
  xpos = 0.0;
  ypos = 0.0;
  ypos  = fgkClength[ilayer][0] + fgkClength[ilayer][1] + fgkClength[ilayer][2]/2.0;
  for (Int_t ic = 0; ic < istack; ic++) {
    ypos -= fgkClength[ilayer][ic];
  }
  ypos -= fgkClength[ilayer][istack]/2.0;
  zpos  = fgkVrocsm + fgkSMpltT + fgkCHsv/2.0 - fgkSheight/2.0
        + ilayer * (fgkCH + fgkVspace);
  TGeoVolume *sm1 = gGeoManager->GetVolume("UTI1");
  TGeoVolume *sm2 = gGeoManager->GetVolume("UTI2");
  TGeoVolume *sm3 = gGeoManager->GetVolume("UTI3");
  TGeoVolume *sm4 = gGeoManager->GetVolume("UTI4");
  sm1->AddNode(roc,1,new TGeoTranslation(xpos,ypos,zpos));
  sm2->AddNode(roc,1,new TGeoTranslation(xpos,ypos,zpos));
  if (istack != 2) {
    // w/o middle stack
    sm3->AddNode(roc,1,new TGeoTranslation(xpos,ypos,zpos));
  }
  if (!((ilayer == 4) && (istack == 4))) {
    // Sector 17 w/o L4S4 chamber
    sm4->AddNode(roc,1,new TGeoTranslation(xpos,ypos,zpos));
  }

} 

//_____________________________________________________________________________
Bool_t AliTRDgeometry::RotateBack(Int_t det
                                , const Double_t * const loc
                                , Double_t *glb) const
{
  //
  // Rotates a chambers to transform the corresponding local frame
  // coordinates <loc> into the coordinates of the ALICE restframe <glb>.
  //

  Int_t   sector = GetSector(det);
  Float_t phi = 2.0 * TMath::Pi() /  (Float_t) fgkNsector * ((Float_t) sector + 0.5);

  glb[0] = loc[0] * TMath::Cos(phi) - loc[1] * TMath::Sin(phi);
  glb[1] = loc[0] * TMath::Sin(phi) + loc[1] * TMath::Cos(phi);
  glb[2] = loc[2];

  return kTRUE;

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetDetectorSec(Int_t layer, Int_t stack)
{
  //
  // Convert plane / stack into detector number for one single sector
  //

  return (layer + stack * fgkNlayer);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetDetector(Int_t layer, Int_t stack, Int_t sector)
{
  //
  // Convert layer / stack / sector into detector number
  //

  return (layer + stack * fgkNlayer + sector * fgkNlayer * fgkNstack);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetLayer(Int_t det)
{
  //
  // Reconstruct the layer number from the detector number
  //

  return ((Int_t) (det % fgkNlayer));

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetStack(Int_t det)
{
  //
  // Reconstruct the stack number from the detector number
  //

  return ((Int_t) (det % (fgkNlayer * fgkNstack)) / fgkNlayer);

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetStack(Double_t z, Int_t layer)
{
  //
  // Reconstruct the chamber number from the z position and layer number
  //
  // The return function has to be protected for positiveness !!
  //

  if ((layer <          0) || 
      (layer >= fgkNlayer)) return -1;
	
  Int_t    istck = fgkNstack;
  Double_t zmin  = 0.0;
  Double_t zmax  = 0.0;

  do {
    istck--;
    if (istck < 0) break;
    AliTRDpadPlane *pp = GetPadPlane(layer,istck);
    zmax  = pp->GetRow0();
    Int_t nrows = pp->GetNrows();
    zmin = zmax -         2 * pp->GetLengthOPad() 
                - (nrows-2) * pp->GetLengthIPad() 
                - (nrows-1) * pp->GetRowSpacing();
  } while((z < zmin) || (z > zmax));
  
  return istck;

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetSector(Int_t det)
{
  //
  // Reconstruct the sector number from the detector number
  //

  return ((Int_t) (det / (fgkNlayer * fgkNstack)));

}

//_____________________________________________________________________________
AliTRDpadPlane *AliTRDgeometry::GetPadPlane(Int_t layer, Int_t stack)
{
  //
  // Returns the pad plane for a given plane <pl> and stack <st> number
  //

  if (!fgPadPlaneArray) {
    CreatePadPlaneArray();
  }

  Int_t ipp = GetDetectorSec(layer,stack);
  return ((AliTRDpadPlane *) fgPadPlaneArray->At(ipp));

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetRowMax(Int_t layer, Int_t stack, Int_t /*sector*/)
{
  //
  // Returns the number of rows on the pad plane
  //

  return GetPadPlane(layer,stack)->GetNrows();

}

//_____________________________________________________________________________
Int_t AliTRDgeometry::GetColMax(Int_t layer)
{
  //
  // Returns the number of rows on the pad plane
  //

  return GetPadPlane(layer,0)->GetNcols();

}

//_____________________________________________________________________________
Double_t AliTRDgeometry::GetRow0(Int_t layer, Int_t stack, Int_t /*sector*/)
{
  //
  // Returns the position of the border of the first pad in a row
  //

  return GetPadPlane(layer,stack)->GetRow0();

}

//_____________________________________________________________________________
Double_t AliTRDgeometry::GetCol0(Int_t layer)
{
  //
  // Returns the position of the border of the first pad in a column
  //

  return GetPadPlane(layer,0)->GetCol0();

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::CreateClusterMatrixArray()
{
  //
  // Create the matrices to transform cluster coordinates from the 
  // local chamber system to the tracking coordinate system
  //

  if (!gGeoManager) {
    return kFALSE;
  }

  if(fgClusterMatrixArray)
    return kTRUE;

  TString volPath;
  TString vpStr   = "ALIC_1/B077_1/BSEGMO";
  TString vpApp1  = "_1/BTRD";
  TString vpApp2  = "_1";
  TString vpApp3a = "/UTR1_1/UTS1_1/UTI1_1";
  TString vpApp3b = "/UTR2_1/UTS2_1/UTI2_1";
  TString vpApp3c = "/UTR3_1/UTS3_1/UTI3_1";
  TString vpApp3d = "/UTR4_1/UTS4_1/UTI4_1";

  fgClusterMatrixArray = new TObjArray(kNdet);
  AliAlignObjParams o;

  for (Int_t iLayer = AliGeomManager::kTRD1; iLayer <= AliGeomManager::kTRD6; iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer); iModule++) {
      
      Int_t        isector   = iModule/Nstack();
      Int_t        istack    = iModule%Nstack();
      Int_t        iLayerTRD = iLayer - AliGeomManager::kTRD1;
      Int_t        lid       = GetDetector(iLayerTRD,istack,isector);    

      // Check for disabled supermodules
      volPath  = vpStr;
      volPath += isector;
      volPath += vpApp1;
      volPath += isector;
      volPath += vpApp2;
      switch (isector) {
      case 17:
        if ((istack == 4) && (iLayerTRD == 4)) {
	  continue;
	}
        volPath += vpApp3d;
        break;
      case 13:
      case 14:
      case 15:
        // Check for holes in from of PHOS
        if (istack == 2) {
          continue;
	}
        volPath += vpApp3c;
        break;
      case 11:
      case 12:
        volPath += vpApp3b;
        break;
      default:
        volPath += vpApp3a;
      };
      if (!gGeoManager->CheckPath(volPath)) {
        //AliInfo(Form("Path not found in geometry: %s",volPath.Data()));
	continue;
      }

      UShort_t     volid   = AliGeomManager::LayerToVolUID(iLayer,iModule);
      const char  *symname = AliGeomManager::SymName(volid);
      TGeoPNEntry *pne     = gGeoManager->GetAlignableEntry(symname);
      const char  *path    = symname;
      if (pne) {
        path = pne->GetTitle();
      }
      else {
	continue;
      }
      if (!strstr(path,"ALIC")) {
        AliDebugClass(1,Form("Not a valid path: %s\n",path));
        continue;
      }
      if (!gGeoManager->cd(path)) {
        AliErrorClass(Form("Cannot go to path: %s\n",path));
        continue;
      }
      TGeoHMatrix *m = gGeoManager->GetCurrentMatrix();
      
      TGeoRotation mchange; 
      mchange.RotateY(90); 
      mchange.RotateX(90);

      //
      // Cluster transformation matrix
      //
      TGeoHMatrix  rotMatrix(mchange.Inverse());
      rotMatrix.MultiplyLeft(m);
      Double_t sectorAngle = 20.0 * (isector % 18) + 10.0;
      TGeoHMatrix  rotSector;
      rotSector.RotateZ(sectorAngle);
      rotMatrix.MultiplyLeft(&rotSector.Inverse());

      fgClusterMatrixArray->AddAt(new TGeoHMatrix(rotMatrix),lid);       

    }    
  }

  return kTRUE;

}

//_____________________________________________________________________________
TGeoHMatrix *AliTRDgeometry::GetClusterMatrix(Int_t det)
{
  //
  // Returns the cluster transformation matrix for a given detector
  //

  if (!fgClusterMatrixArray) {
    if (!CreateClusterMatrixArray()) {
      return NULL;
    }
  }  
  return (TGeoHMatrix *) fgClusterMatrixArray->At(det);

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::ChamberInGeometry(Int_t det)
{
  //
  // Checks whether the given detector is part of the current geometry
  //

  if (!GetClusterMatrix(det)) {
    return kFALSE;
  }
  else {
    return kTRUE;
  }

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::IsHole(Int_t /*la*/, Int_t st, Int_t se) const
{
  //
  // Checks for holes in front of PHOS
  //

  if (((se == 13) || (se == 14) || (se == 15)) && 
      (st == 2)) {
    return kTRUE; 
  }

  return kFALSE;

}

//_____________________________________________________________________________
Bool_t AliTRDgeometry::IsOnBoundary(Int_t det, Float_t y, Float_t z, Float_t eps) const
{
  //
  // Checks whether position is at the boundary of the sensitive volume 
  //

  Int_t ly = GetLayer(det);
  if ((ly <          0) || 
      (ly >= fgkNlayer)) return kTRUE;
	
  Int_t stk = GetStack(det);
  if ((stk <          0) || 
      (stk >= fgkNstack)) return kTRUE;

  AliTRDpadPlane *pp = (AliTRDpadPlane*) fgPadPlaneArray->At(GetDetectorSec(ly, stk));
  if(!pp) return kTRUE;

  Double_t max  = pp->GetRow0();
  Int_t n = pp->GetNrows();
  Double_t min = max - 2 * pp->GetLengthOPad() 
                 - (n-2) * pp->GetLengthIPad() 
                 - (n-1) * pp->GetRowSpacing();
  if(z < min+eps || z > max-eps){ 
    //printf("z : min[%7.2f (%7.2f)] %7.2f max[(%7.2f) %7.2f]\n", min, min+eps, z, max-eps, max);
    return kTRUE;
  }
  min  = pp->GetCol0();
  n = pp->GetNcols();
  max = min +2 * pp->GetWidthOPad() 
       + (n-2) * pp->GetWidthIPad() 
       + (n-1) * pp->GetColSpacing();
  if(y < min+eps || y > max-eps){ 
    //printf("y : min[%7.2f (%7.2f)] %7.2f max[(%7.2f) %7.2f]\n", min, min+eps, y, max-eps, max);
    return kTRUE;
  }

  return kFALSE;

}
