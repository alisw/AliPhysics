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
Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.2  2000/05/08 14:46:44  cblume
Include options SetPHOShole() and SetRICHhole()

Revision 1.1.4.1  2000/04/27 12:46:04  cblume
Corrected bug in full geometry

Revision 1.1  2000/02/28 19:01:15  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD geometry for the spaceframe without holes                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDgeometryFull.h"

ClassImp(AliTRDgeometryFull)

//_____________________________________________________________________________
AliTRDgeometryFull::AliTRDgeometryFull():AliTRDgeometry()
{
  //
  // AliTRDgeometryFull default constructor
  //

  Init();

}

//_____________________________________________________________________________
AliTRDgeometryFull::~AliTRDgeometryFull()
{
  //
  // AliTRDgeometryFull destructor
  //

}

//_____________________________________________________________________________
void AliTRDgeometryFull::Init()
{
  //
  // Initializes the geometry parameter
  //

  Int_t iplan;

  fPHOShole = kFALSE;
  fRICHhole = kFALSE;

  // The length of the inner chambers
  for (iplan = 0; iplan < kNplan; iplan++) 
    fClengthI[iplan] = 110.0;

  // The length of the middle chambers
  fClengthM1[0] = 123.5;
  fClengthM1[1] = 131.0;
  fClengthM1[2] = 138.5;
  fClengthM1[3] = 146.0;
  fClengthM1[4] = 153.0;
  fClengthM1[5] = 160.5;

  fClengthM2[0] = 123.5 - 7.0;
  fClengthM2[1] = 131.0 - 7.0;
  fClengthM2[2] = 138.5 - 7.0;
  fClengthM2[3] = 146.0 - 7.0;
  fClengthM2[4] = 153.0 - 7.0;
  fClengthM2[5] = 160.4 - 7.0;

  // The length of the outer chambers
  fClengthO1[0] = 123.5;
  fClengthO1[1] = 131.0;
  fClengthO1[2] = 134.5;
  fClengthO1[3] = 142.0;
  fClengthO1[4] = 142.0;
  fClengthO1[5] = 134.5;

  fClengthO2[0] = 123.5;
  fClengthO2[1] = 131.0;
  fClengthO2[2] = 134.5;
  fClengthO2[3] = 142.0;
  fClengthO2[4] = 142.0;
  fClengthO2[5] = 134.5;

  fClengthO3[0] =  86.5;
  fClengthO3[1] = 101.5;
  fClengthO3[2] = 112.5;
  fClengthO3[3] = 127.5;
  fClengthO3[4] = 134.5;
  fClengthO3[5] = 134.5;

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

  // The pad row (z-direction)
  for (iplan = 0; iplan < kNplan; iplan++) {

    for (Int_t isect = 0; isect < kNsect; isect++) {
      Float_t clengthI = fClengthI[iplan];
      Float_t clengthM = fClengthM1[iplan];
      Float_t clengthO = fClengthO1[iplan];
      switch (isect) {
      case 12:
      case 13:
      case 14:
      case 15:
      case 16:
        clengthM = fClengthM2[iplan];
        clengthO = fClengthO2[iplan];
        break;
      case 4:
      case 5:
      case 6:
        clengthO = fClengthO3[iplan];
        break;
      };
      fRowMax[iplan][0][isect] = 1 + TMath::Nint((clengthO - 2. * kCcthick) 
                                                           / fRowPadSize - 0.5);
      fRowMax[iplan][1][isect] = 1 + TMath::Nint((clengthM - 2. * kCcthick) 
                                                           / fRowPadSize - 0.5);
      fRowMax[iplan][2][isect] = 1 + TMath::Nint((clengthI - 2. * kCcthick) 
                                                           / fRowPadSize - 0.5);
      fRowMax[iplan][3][isect] = 1 + TMath::Nint((clengthM - 2. * kCcthick) 
                                                           / fRowPadSize - 0.5);
      fRowMax[iplan][4][isect] = 1 + TMath::Nint((clengthO - 2. * kCcthick) 
                                                           / fRowPadSize - 0.5);
      fRow0[iplan][0][isect]   = -clengthI/2. - clengthM - clengthO + kCcthick; 
      fRow0[iplan][1][isect]   = -clengthI/2. - clengthM            + kCcthick;
      fRow0[iplan][2][isect]   = -clengthI/2.                       + kCcthick;
      fRow0[iplan][3][isect]   =  clengthI/2.                       + kCcthick; 
      fRow0[iplan][4][isect]   =  clengthI/2. + clengthM            + kCcthick; 
    }

  }

}

//_____________________________________________________________________________
void AliTRDgeometryFull::CreateGeometry(Int_t *idtmed)
{
  //
  // Create the TRD geometry without hole
  //

  Int_t iplan;

  const Int_t kNparTrd = 4;
  const Int_t kNparCha = 3;

  Float_t parTrd[kNparTrd];
  Float_t parCha[kNparCha];

  Float_t xpos, ypos, zpos;

  AliTRDgeometry::CreateGeometry(idtmed);

  // The TRD mother volume for one sector (Air), full length in z-direction
  parTrd[0] = kSwidth1/2.;
  parTrd[1] = kSwidth2/2.;
  parTrd[2] = kSlenTR1/2.;
  parTrd[3] = kSheight/2.;
  gMC->Gsvolu("TRD1","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  
  // The TRD mother volume for one sector (Air), leaving hole for PHOS
  if (fPHOShole) {
    gMC->Gsvolu("TRD2","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  }

  // The TRD mother volume for one sector (Air), leaving hole for RICH
  if (fRICHhole) {
    gMC->Gsvolu("TRD3","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  }  

  // Position the chambers in the TRD mother volume
  for (iplan = 1; iplan <= kNplan; iplan++) {

    Float_t y1 = fClengthM1[iplan-1] - fClengthM2[iplan-1];
    Float_t y2 = fClengthO1[iplan-1] - fClengthO3[iplan-1];

    // The inner chambers ---------------------------------------------------------------

    // the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthI[iplan-1]/2.;
    parCha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAFI",iplan       ,"TRD1",xpos,ypos,zpos,0,"MANY",parCha,kNparCha);

    // the inner part of the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.   - kCathick;
    parCha[1] = fClengthI[iplan-1]/2. - kCathick;
    parCha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAII",iplan       ,"TRD1",xpos,ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthI[iplan-1]/2.;
    parCha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCFI",iplan       ,"TRD1",xpos,ypos,zpos,0,"MANY",parCha,kNparCha);

    // the inner part of the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.   - kCcthick;
    parCha[1] = fClengthI[iplan-1]/2. - kCcthick;
    parCha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCII",iplan       ,"TRD1",xpos,ypos,zpos,0,"ONLY",parCha,kNparCha);

    // The middle chambers --------------------------------------------------------------

    // the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthM1[iplan-1]/2.;
    parCha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAFM",iplan         ,"TRD1",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
    gMC->Gsposp("UAFM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    if (fPHOShole) {
      parCha[0] = fCwidth[iplan-1]/2.;
      parCha[1] = fClengthM2[iplan-1]/2.;
      parCha[2] = kCaframe/2.;
      xpos       = 0.;
      ypos       = fClengthI[iplan-1]/2. + fClengthM2[iplan-1]/2. + y1;
      zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UAFM",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
      gMC->Gsposp("UAFM",iplan+3*kNplan,"TRD2",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    }

    // the inner part of the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.   - kCathick;
    parCha[1] = fClengthM1[iplan-1]/2. - kCathick;
    parCha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAIM",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UAIM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    if (fPHOShole) {
      parCha[0] = fCwidth[iplan-1]/2.    - kCathick;
      parCha[1] = fClengthM2[iplan-1]/2. - kCathick;
      parCha[2] = kCaframe/2.;
      xpos       = 0.;
      ypos       = fClengthI[iplan-1]/2. + fClengthM2[iplan-1]/2. + y1;
      zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UAIM",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
      gMC->Gsposp("UAIM",iplan+3*kNplan,"TRD2",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    }

    // the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthM1[iplan-1]/2.;
    parCha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos       = kCcframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCFM",iplan         ,"TRD1",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
    gMC->Gsposp("UCFM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    if (fPHOShole) {
      parCha[0] = fCwidth[iplan-1]/2.;
      parCha[1] = fClengthM2[iplan-1]/2.;
      parCha[2] = kCcframe/2.;
      xpos       = 0.;
      ypos       = fClengthI[iplan-1]/2. + fClengthM2[iplan-1]/2. + y1;
      zpos       = kCcframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UCFM",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
      gMC->Gsposp("UCFM",iplan+3*kNplan,"TRD2",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    }

    // the inner part of the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.   - kCcthick;
    parCha[1] = fClengthM1[iplan-1]/2. - kCcthick;
    parCha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos       = kCcframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCIM",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UCIM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    if (fPHOShole) {
      parCha[0] = fCwidth[iplan-1]/2.    - kCcthick;
      parCha[1] = fClengthM2[iplan-1]/2. - kCcthick;
      parCha[2] = kCcframe/2.;
      xpos       = 0.;
      ypos       = fClengthI[iplan-1]/2. + fClengthM2[iplan-1]/2. + y1;
      zpos       = kCcframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UCIM",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
      gMC->Gsposp("UCIM",iplan+3*kNplan,"TRD2",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    }

    // The outer chambers ---------------------------------------------------------------

    // the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthO1[iplan-1]/2.;
    parCha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO1[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAFO",iplan         ,"TRD1",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
    gMC->Gsposp("UAFO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    if (fPHOShole) {
      parCha[0] = fCwidth[iplan-1]/2.;
      parCha[1] = fClengthO2[iplan-1]/2.;
      parCha[2] = kCaframe/2.;
      xpos       = 0.;
      ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO2[iplan-1]/2.;
      zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UAFO",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
      gMC->Gsposp("UAFO",iplan+3*kNplan,"TRD2",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    }
    if (fRICHhole) {
      parCha[0] = fCwidth[iplan-1]/2.;
      parCha[1] = fClengthO3[iplan-1]/2.;
      parCha[2] = kCaframe/2.;
      xpos       = 0.;
      ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO3[iplan-1]/2. + y2;
      zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UAFO",iplan+4*kNplan,"TRD3",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
      gMC->Gsposp("UAFO",iplan+5*kNplan,"TRD3",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    }

    // the inner part of the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.   - kCathick;
    parCha[1] = fClengthO1[iplan-1]/2. - kCathick;
    parCha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO1[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAIO",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UAIO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    if (fPHOShole) {
      parCha[0] = fCwidth[iplan-1]/2.    - kCathick;
      parCha[1] = fClengthO2[iplan-1]/2. - kCathick;
      parCha[2] = kCaframe/2.;
      xpos       = 0.;
      ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO2[iplan-1]/2.;
      zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UAIO",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
      gMC->Gsposp("UAIO",iplan+3*kNplan,"TRD2",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    }
    if (fRICHhole) {
      parCha[0] = fCwidth[iplan-1]/2.    - kCathick;
      parCha[1] = fClengthO3[iplan-1]/2. - kCathick;
      parCha[2] = kCaframe/2.;
      xpos       = 0.;
      ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO3[iplan-1]/2. + y2;
      zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UAIO",iplan+4*kNplan,"TRD3",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
      gMC->Gsposp("UAIO",iplan+5*kNplan,"TRD3",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    }

    // the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthO1[iplan-1]/2.;
    parCha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO1[iplan-1]/2.;
    zpos       = kCcframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCFO",iplan,         "TRD1",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
    gMC->Gsposp("UCFO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    if (fPHOShole) {
      parCha[0] = fCwidth[iplan-1]/2.;
      parCha[1] = fClengthO2[iplan-1]/2.;
      parCha[2] = kCcframe/2.;
      xpos       = 0.;
      ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO2[iplan-1]/2.;
      zpos       = kCcframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UCFO",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
      gMC->Gsposp("UCFO",iplan+3*kNplan,"TRD2",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    }
    if (fRICHhole) {
      parCha[0] = fCwidth[iplan-1]/2.;
      parCha[1] = fClengthO3[iplan-1]/2.;
      parCha[2] = kCcframe/2.;
      xpos       = 0.;
      ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO3[iplan-1]/2. + y2;
      zpos       = kCcframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UCFO",iplan+4*kNplan,"TRD3",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
      gMC->Gsposp("UCFO",iplan+5*kNplan,"TRD3",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    }

    // the inner part of the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.   - kCcthick;
    parCha[1] = fClengthO1[iplan-1]/2. - kCcthick;
    parCha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO1[iplan-1]/2.;
    zpos       = kCcframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCIO",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UCIO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    if (fPHOShole) {
      parCha[0] = fCwidth[iplan-1]/2.    - kCcthick;
      parCha[1] = fClengthO2[iplan-1]/2. - kCcthick;
      parCha[2] = kCcframe/2.;
      xpos       = 0.;
      ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO2[iplan-1]/2.;
      zpos       = kCcframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UCIO",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
      gMC->Gsposp("UCIO",iplan+3*kNplan,"TRD2",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    }
    if (fRICHhole) {
      parCha[0] = fCwidth[iplan-1]/2.    - kCcthick;
      parCha[1] = fClengthO3[iplan-1]/2. - kCcthick;
      parCha[2] = kCcframe/2.;
      xpos       = 0.;
      ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO3[iplan-1]/2. + y2;
      zpos       = kCcframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UCIO",iplan+4*kNplan,"TRD3",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
      gMC->Gsposp("UCIO",iplan+5*kNplan,"TRD3",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    }

  }

  xpos     = 0.;
  ypos     = 0.;
  zpos     = 0.;
  gMC->Gspos("TRD1",1,"BTR1",xpos,ypos,zpos,0,"ONLY");
  if (fPHOShole) 
    gMC->Gspos("TRD2",2,"BTR2",xpos,ypos,zpos,0,"ONLY");
  else
    gMC->Gspos("TRD1",2,"BTR2",xpos,ypos,zpos,0,"ONLY");
  if (fRICHhole)
    gMC->Gspos("TRD3",3,"BTR3",xpos,ypos,zpos,0,"ONLY");
  else
    gMC->Gspos("TRD1",3,"BTR3",xpos,ypos,zpos,0,"ONLY");

}
