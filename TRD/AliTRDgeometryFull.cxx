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
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD geometry without holes                                               //
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

}

//_____________________________________________________________________________
void AliTRDgeometryFull::Init()
{
  //
  // Initializes the geometry parameter
  //

  Int_t iplan;

  // The length of the inner chambers
  for (iplan = 0; iplan < kNplan; iplan++) 
    fClengthI[iplan] = 110.0;

  // The length of the middle chambers
  fClengthM[0] = 123.5;
  fClengthM[1] = 131.0;
  fClengthM[2] = 138.5;
  fClengthM[3] = 146.0;
  fClengthM[4] = 153.0;
  fClengthM[5] = 160.5;

  // The length of the outer chambers
  fClengthO[0] = 123.5;
  fClengthO[1] = 131.0;
  fClengthO[2] = 134.5;
  fClengthO[3] = 142.0;
  fClengthO[4] = 142.0;
  fClengthO[5] = 134.5;

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
      Float_t clengthM = fClengthM[iplan];
      Float_t clengthO = fClengthO[iplan];
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

  const Int_t npar_trd = 4;
  const Int_t npar_cha = 3;

  Float_t par_trd[npar_trd];
  Float_t par_cha[npar_cha];

  Float_t xpos, ypos, zpos;

  AliTRDgeometry::CreateGeometry(idtmed);

  // The TRD mother volume for one sector (Air) (dimensions identical to BTR1)
  par_trd[0] = kSwidth1/2.;
  par_trd[1] = kSwidth2/2.;
  par_trd[2] = kSlenTR1/2.;
  par_trd[3] = kSheight/2.;
  gMC->Gsvolu("TRD1","TRD1",idtmed[1302-1],par_trd,npar_trd);
  
  // Position the chambers in the TRD mother volume
  for (iplan = 1; iplan <= kNplan; iplan++) {

    // The inner chambers ---------------------------------------------------------------

    // the aluminum frame
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthI[iplan-1]/2.;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAFI",iplan       ,"TRD1",xpos,ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the aluminum frame
    par_cha[0] = fCwidth[iplan-1]/2.   - kCathick;
    par_cha[1] = fClengthI[iplan-1]/2. - kCathick;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAII",iplan       ,"TRD1",xpos,ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // the carbon frame
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthI[iplan-1]/2.;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCFI",iplan       ,"TRD1",xpos,ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the carbon frame
    par_cha[0] = fCwidth[iplan-1]/2.   - kCcthick;
    par_cha[1] = fClengthI[iplan-1]/2. - kCcthick;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCII",iplan       ,"TRD1",xpos,ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // The middle chambers --------------------------------------------------------------

    // the aluminum frame
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthM[iplan-1]/2.;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2.  + fClengthM[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAFM",iplan         ,"TRD1",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);
    gMC->Gsposp("UAFM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the aluminum frame
    par_cha[0] = fCwidth[iplan-1]/2.    - kCathick;
    par_cha[1] = fClengthM[iplan-1]/2. - kCathick;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2.  + fClengthM[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAIM",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
    gMC->Gsposp("UAIM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // the carbon frame
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthM[iplan-1]/2.;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]/2.;
    zpos       = kCcframe/2.           - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCFM",iplan         ,"TRD1",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);
    gMC->Gsposp("UCFM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the carbon frame
    par_cha[0] = fCwidth[iplan-1]/2.    - kCcthick;
    par_cha[1] = fClengthM[iplan-1]/2. - kCcthick;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]/2.;
    zpos       = kCcframe/2.           - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCIM",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
    gMC->Gsposp("UCIM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // The outer chambers ---------------------------------------------------------------

    // the aluminum frame
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthO[iplan-1]/2.;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]    + fClengthO[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAFO",iplan         ,"TRD1",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);
    gMC->Gsposp("UAFO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the aluminum frame
    par_cha[0] = fCwidth[iplan-1]/2.    - kCathick;
    par_cha[1] = fClengthO[iplan-1]/2. - kCathick;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]    + fClengthO[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAIO",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
    gMC->Gsposp("UAIO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // the carbon frame
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthO[iplan-1]/2.;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]    + fClengthO[iplan-1]/2.;
    zpos       = kCcframe/2.           - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCFO",iplan,         "TRD1",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);
    gMC->Gsposp("UCFO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the carbon frame
    par_cha[0] = fCwidth[iplan-1]/2.    - kCcthick;
    par_cha[1] = fClengthO[iplan-1]/2. - kCcthick;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]    + fClengthO[iplan-1]/2.;
    zpos       = kCcframe/2.           - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCIO",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
    gMC->Gsposp("UCIO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",par_cha,npar_cha);

  }

  xpos     = 0.;
  ypos     = 0.;
  zpos     = 0.;
  gMC->Gspos("TRD1",1,"BTR1",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("TRD1",2,"BTR2",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("TRD1",3,"BTR3",xpos,ypos,zpos,0,"ONLY");

}
