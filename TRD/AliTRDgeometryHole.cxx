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
Revision 1.4  2000/11/01 14:53:21  cblume
Merge with TRD-develop

Revision 1.1.4.4  2000/10/15 23:40:01  cblume
Remove AliTRDconst

Revision 1.1.4.3  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.4.2  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.1.4.1  2000/09/22 14:43:41  cblume
Allow the pad/timebin-dimensions to be changed after initialization

Revision 1.3  2000/10/02 21:28:19  fca
Removal of useless dependecies via forward declarations

Revision 1.2  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.1  2000/02/28 19:01:42  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD geometry with holes                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliMC.h"

#include "AliTRDgeometryHole.h"

ClassImp(AliTRDgeometryHole)

//_____________________________________________________________________________
AliTRDgeometryHole::AliTRDgeometryHole():AliTRDgeometry()
{
  //
  // AliTRDgeometryHole default constructor
  //

  Init();

}

//_____________________________________________________________________________
AliTRDgeometryHole::~AliTRDgeometryHole()
{
  //
  // AliTRDgeometryHole destructor
  //

}

//_____________________________________________________________________________
void AliTRDgeometryHole::Init()
{
  //
  // Initializes the geometry parameter
  //

  Int_t iplan;

  // The length of the inner chambers
  for (iplan = 0; iplan < fgkNplan; iplan++) 
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
  SetNRowPad();

}

//_____________________________________________________________________________
void AliTRDgeometryHole::SetNRowPad(Int_t p, Int_t c, Int_t npad)
{
  //
  // Redefines the number of pads in raw direction for
  // a given plane and chamber number
  //

  Float_t clengthI = fClengthI[p];
  Float_t clengthM = fClengthM1[p];
  Float_t clengthO = fClengthO1[p];
  
  for (Int_t iSect = 0; iSect < fgkNsect; iSect++) {

    fRowMax[p][c][iSect] = npad;
    if (c == 2) {
      fRowPadSize[p][c][iSect] = (clengthI - 2. * fgkCcthick)
                               / fRowMax[p][c][iSect];
    }
    if ((c == 1) || (c == 3)) {
      fRowPadSize[p][c][iSect] = (clengthM - 2. * fgkCcthick)
                               / fRowMax[p][c][iSect];
    }
    if ((c == 0) || (c == 4)) {
      fRowPadSize[p][c][iSect] = (clengthO - 2. * fgkCcthick)
                               / fRowMax[p][c][iSect];
    }

  }

}

//_____________________________________________________________________________
void AliTRDgeometryHole::SetNRowPad()
{
  //
  // Defines the pad size in row direction
  //

  for (Int_t iPlan = 0; iPlan < fgkNplan; iPlan++) {

    Float_t clengthI = fClengthI[iPlan];
    Float_t clengthM = fClengthM1[iPlan];
    Float_t clengthO = fClengthO1[iPlan];

    for (Int_t iSect = 0; iSect < fgkNsect; iSect++) {

      fRow0[iPlan][0][iSect]   = -clengthI/2. - clengthM - clengthO + fgkCcthick; 
      fRow0[iPlan][1][iSect]   = -clengthI/2. - clengthM            + fgkCcthick;
      fRow0[iPlan][2][iSect]   = -clengthI/2.                       + fgkCcthick;
      fRow0[iPlan][3][iSect]   =  clengthI/2.                       + fgkCcthick; 
      fRow0[iPlan][4][iSect]   =  clengthI/2. + clengthM            + fgkCcthick; 

      for (Int_t iCham = 0; iCham < fgkNcham; iCham++) {

        if (iCham == 2) {
          fRowMax[iPlan][iCham][iSect] = 18;
          fRowPadSize[iPlan][iCham][iSect] = (clengthI - 2. * fgkCcthick)
                                           / fRowMax[iPlan][iCham][iSect];
        }
        if ((iCham == 1) || (iCham == 3)) {
          fRowMax[iPlan][iCham][iSect] = 24;
          fRowPadSize[iPlan][iCham][iSect] = (clengthM - 2. * fgkCcthick)
                                           / fRowMax[iPlan][iCham][iSect];
        }
        if ((iCham == 0) || (iCham == 4)) {
          if (iPlan <  4) {
            fRowMax[iPlan][iCham][iSect] = 24;
          }
          if (iPlan == 4) {
            fRowMax[iPlan][iCham][iSect] = 22;
          }
          if (iPlan == 5) {
            fRowMax[iPlan][iCham][iSect] = 20;
          }
          fRowPadSize[iPlan][iCham][iSect] = (clengthO - 2. * fgkCcthick)
                                           / fRowMax[iPlan][iCham][iSect];
        }

      }
    }
  }

}

//_____________________________________________________________________________
void AliTRDgeometryHole::CreateGeometry(Int_t *idtmed)
{
  //
  // Create the TRD geometry with hole
  //

  Int_t iplan;

  const Int_t kNparTrd = 4;
  const Int_t kNparCha = 3;
  const Int_t kNplan   = fgkNplan;

  Float_t parTrd[kNparTrd];
  Float_t parCha[kNparCha];

  Float_t xpos, ypos, zpos;

  AliTRDgeometry::CreateGeometry(idtmed);

  // The TRD mother volume for one sector (Air) (dimensions identical to BTR1)
  parTrd[0] = fgkSwidth1/2.;
  parTrd[1] = fgkSwidth2/2.;
  parTrd[2] = fgkSlenTR1/2.;
  parTrd[3] = fgkSheight/2.;
  gMC->Gsvolu("TRD1","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  
  // The TRD mother volume for one sector (Air) (dimensions identical to BTR2)
  parTrd[0] = fgkSwidth1/2.;
  parTrd[1] = fgkSwidth2/2.;
  parTrd[2] = fgkSlenTR2/2.;
  parTrd[3] = fgkSheight/2.;
  gMC->Gsvolu("TRD2","TRD1",idtmed[1302-1],parTrd,kNparTrd);

  // The TRD mother volume for one sector (Air) (dimensions identical to BTR3)
  parTrd[0] = fgkSwidth1/2.;
  parTrd[1] = fgkSwidth2/2.;
  parTrd[2] = fgkSlenTR3/2.;
  parTrd[3] = fgkSheight/2.;
  gMC->Gsvolu("TRD3","TRD1",idtmed[1302-1],parTrd,kNparTrd);

  // Position the chambers in the TRD mother volume
  for (iplan = 1; iplan <= kNplan; iplan++) {

    // The inner chambers ---------------------------------------------------------------

    // the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthI[iplan-1]/2.;
    parCha[2] = fgkCaframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = fgkCheight - fgkCaframe/2. - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAFI",iplan       ,"TRD1",xpos,ypos,zpos,0,"MANY",parCha,kNparCha);

    // the inner part of the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.   - fgkCathick;
    parCha[1] = fClengthI[iplan-1]/2. - fgkCathick;
    parCha[2] = fgkCaframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = fgkCheight - fgkCaframe/2. - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAII",iplan       ,"TRD1",xpos,ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthI[iplan-1]/2.;
    parCha[2] = fgkCcframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = fgkCcframe/2.            - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCFI",iplan       ,"TRD1",xpos,ypos,zpos,0,"MANY",parCha,kNparCha);

    // the inner part of the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.   - fgkCcthick;
    parCha[1] = fClengthI[iplan-1]/2. - fgkCcthick;
    parCha[2] = fgkCcframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = fgkCcframe/2.            - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCII",iplan       ,"TRD1",xpos,ypos,zpos,0,"ONLY",parCha,kNparCha);

    // The middle chambers --------------------------------------------------------------

    // the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthM1[iplan-1]/2.;
    parCha[2] = fgkCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2.+ fClengthM1[iplan-1]/2.;
    zpos       = fgkCheight - fgkCaframe/2. - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAFM",iplan         ,"TRD1",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
    gMC->Gsposp("UAFM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthM2[iplan-1]/2.;
    parCha[2] = fgkCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthM2[iplan-1]/2. - fgkSlenTR2/2.;
    zpos       = fgkCheight - fgkCaframe/2. - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAFM",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);

    // the inner part of the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.    - fgkCathick;
    parCha[1] = fClengthM1[iplan-1]/2. - fgkCathick;
    parCha[2] = fgkCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos       = fgkCheight - fgkCaframe/2. - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAIM",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UAIM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    parCha[0] = fCwidth[iplan-1]/2.    - fgkCathick;
    parCha[1] = fClengthM2[iplan-1]/2. - fgkCathick;
    parCha[2] = fgkCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthM2[iplan-1]/2. - fgkSlenTR2/2.;
    zpos       = fgkCheight - fgkCaframe/2. - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAIM",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthM1[iplan-1]/2.;
    parCha[2] = fgkCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos       = fgkCcframe/2.         - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCFM",iplan         ,"TRD1",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
    gMC->Gsposp("UCFM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthM2[iplan-1]/2.;
    parCha[2] = fgkCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthM2[iplan-1]/2. - fgkSlenTR2/2.;
    zpos       = fgkCcframe/2.          - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCFM",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);

    // the inner part of the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.    - fgkCcthick;
    parCha[1] = fClengthM1[iplan-1]/2. - fgkCcthick;
    parCha[2] = fgkCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos       = fgkCcframe/2.         - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCIM",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UCIM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    parCha[0] = fCwidth[iplan-1]/2.    - fgkCcthick;
    parCha[1] = fClengthM2[iplan-1]/2. - fgkCcthick;
    parCha[2] = fgkCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthM2[iplan-1]/2. - fgkSlenTR2/2.;
    zpos       = fgkCcframe/2.          - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCIM",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);

    // The outer chambers ---------------------------------------------------------------

    // the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthO1[iplan-1]/2.;
    parCha[2] = fgkCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]    + fClengthO1[iplan-1]/2.;
    zpos       = fgkCheight - fgkCaframe/2. - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAFO",iplan         ,"TRD1",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
    gMC->Gsposp("UAFO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthO2[iplan-1]/2.;
    parCha[2] = fgkCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthM2[iplan-1]   + fClengthO2[iplan-1]/2. - fgkSlenTR2/2.;
    zpos       = fgkCheight - fgkCaframe/2. - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAFO",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthO3[iplan-1]/2.;
    parCha[2] = fgkCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthO3[iplan-1]/2. - fgkSlenTR3/2.;
    zpos       = fgkCheight - fgkCaframe/2. - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAFO",iplan+4*kNplan,"TRD3",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);

    // the inner part of the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.    - fgkCathick;
    parCha[1] = fClengthO1[iplan-1]/2. - fgkCathick;
    parCha[2] = fgkCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]    + fClengthO1[iplan-1]/2.;
    zpos       = fgkCheight - fgkCaframe/2. - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAIO",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UAIO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    parCha[0] = fCwidth[iplan-1]/2.    - fgkCathick;
    parCha[1] = fClengthO2[iplan-1]/2. - fgkCathick;
    parCha[2] = fgkCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthM2[iplan-1]   + fClengthO2[iplan-1]/2. - fgkSlenTR2/2.;
    zpos       = fgkCheight - fgkCaframe/2. - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAIO",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    parCha[0] = fCwidth[iplan-1]/2.    - fgkCathick;
    parCha[1] = fClengthO3[iplan-1]/2. - fgkCathick;
    parCha[2] = fgkCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthO3[iplan-1]/2. - fgkSlenTR3/2.;
    zpos       = fgkCheight - fgkCaframe/2. - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAIO",iplan+4*kNplan,"TRD3",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthO1[iplan-1]/2.;
    parCha[2] = fgkCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]    + fClengthO1[iplan-1]/2.;
    zpos       = fgkCcframe/2.         - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCFO",iplan,         "TRD1",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
    gMC->Gsposp("UCFO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",parCha,kNparCha);
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthO2[iplan-1]/2.;
    parCha[2] = fgkCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthM2[iplan-1]    + fClengthO2[iplan-1]/2. - fgkSlenTR2/2.;
    zpos       = fgkCcframe/2.          - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCFO",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthO3[iplan-1]/2.;
    parCha[2] = fgkCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthO3[iplan-1]/2. - fgkSlenTR3/2.;
    zpos       = fgkCcframe/2.          - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCFO",iplan+4*kNplan,"TRD3",xpos, ypos,zpos,0,"MANY",parCha,kNparCha);

    // the inner part of the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.    - fgkCcthick;
    parCha[1] = fClengthO1[iplan-1]/2. - fgkCcthick;
    parCha[2] = fgkCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]    + fClengthO1[iplan-1]/2.;
    zpos       = fgkCcframe/2.         - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCIO",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UCIO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);
    parCha[0] = fCwidth[iplan-1]/2.    - fgkCcthick;
    parCha[1] = fClengthO2[iplan-1]/2. - fgkCcthick;
    parCha[2] = fgkCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthM2[iplan-1]    + fClengthO2[iplan-1]/2. - fgkSlenTR2/2.;
    zpos       = fgkCcframe/2.          - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCIO",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    parCha[0] = fCwidth[iplan-1]/2.    - fgkCcthick;
    parCha[1] = fClengthO3[iplan-1]/2. - fgkCcthick;
    parCha[2] = fgkCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthO3[iplan-1]/2. - fgkSlenTR3/2.;
    zpos       = fgkCcframe/2.          - fgkSheight/2. + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCIO",iplan+4*kNplan,"TRD3",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);

  }

  xpos     = 0.;
  ypos     = 0.;
  zpos     = 0.;
  gMC->Gspos("TRD1",1,"BTR1",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("TRD2",1,"BTR2",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("TRD3",1,"BTR3",xpos,ypos,zpos,0,"ONLY");

}
