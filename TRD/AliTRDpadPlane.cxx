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
//  Describes a pad plane of a TRD ROC                                       //
//                                                                           //
//  Contains the information on pad postions, pad dimensions,                //
//  tilting angle, etc.                                                      //
//  It also provides methods to identify the current pad number from         //
//  global coordinates.                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDpadPlane.h"
#include "AliTRDgeometryFull.h"

ClassImp(AliTRDpadPlane)

//_____________________________________________________________________________
AliTRDpadPlane::AliTRDpadPlane():TObject()
{
  //
  // Default constructor
  //

  fGeo          = 0;

  fPla          = 0;
  fCha          = 0;

  fLength       = 0.0;
  fWidth        = 0.0;
  fLengthRim    = 0.0;
  fWidthRim     = 0.0;
  fLengthOPad   = 0.0;
  fWidthOPad    = 0.0;
  fLengthIPad   = 0.0;
  fWidthIPad    = 0.0;

  fRowSpacing   = 0.0;
  fColSpacing   = 0.0;

  fNrows        = 0;
  fNcols        = 0;

  fPadRow       = 0;
  fPadCol       = 0;

  fTiltingAngle = 0.0;

}

//_____________________________________________________________________________
AliTRDpadPlane::AliTRDpadPlane(Int_t p, Int_t c):TObject()
{
  //
  // Constructor that initializes a given pad plane type
  //

  fGeo = new AliTRDgeometryFull();

  fPla = p;
  fCha = c;

  //fRowSpacing = 0.025;
  //fColSpacing = 0.025;

  fRowSpacing = 0.0;
  fColSpacing = 0.0;

  fLengthRim  = 1.0;
  fWidthRim   = 0.5;

  fNcols      = 144;

  //
  // The pad plane parameter
  //
  switch (p) {
  case 0:
    if (c == 2) {
      // L0C0 type
      fNrows        =  12;
      fLength       = 108.0;
      fWidth        =  94.4;
      fLengthOPad   =   8.0;
      fWidthOPad    =   0.515;
      fLengthIPad   =   9.0;
      fWidthIPad    =   0.635;
      fTiltingAngle =  -2.0;
    }
    else {
      // L0C1 type
      fNrows        =  16;
      fLength       = 122.0;
      fWidth        =  94.4;
      fLengthOPad   =   7.5;
      fWidthOPad    =   0.515;
      fLengthIPad   =   7.5;
      fWidthIPad    =   0.635;
      fTiltingAngle =  -2.0;
    }
    break;
  case 1:
    if (c == 2) {
      // L1C0 type
      fNrows        =  12;
      fLength       = 108.0;
      fWidth        =  96.6;
      fLengthOPad   =   8.0;
      fWidthOPad    =   0.585;
      fLengthIPad   =   9.0;
      fWidthIPad    =   0.665;
      fTiltingAngle =   2.0;
    }
    else {
      // L1C1 type
      fNrows        =  16;
      fLength       = 122.0;
      fWidth        =  96.6;
      fLengthOPad   =   7.5;
      fWidthOPad    =   0.585;
      fLengthIPad   =   7.5;
      fWidthIPad    =   0.665;
      fTiltingAngle =   2.0;
    }
    break;
  case 2:
    if (c == 2) {
      // L2C0 type
      fNrows        =  12;
      fLength       = 108.0;
      fWidth        = 101.1;
      fLengthOPad   =   8.0;
      fWidthOPad    =   0.705;
      fLengthIPad   =   9.0;
      fWidthIPad    =   0.695;
      fTiltingAngle =  -2.0;
    }
    else {
      // L2C1 type
      fNrows        =  16;
      fLength       = 129.0;
      fWidth        = 101.1;
      fLengthOPad   =   7.5;
      fWidthOPad    =   0.705;
      fLengthIPad   =   8.0;
      fWidthIPad    =   0.695;
      fTiltingAngle =  -2.0;
    }
    break;
  case 3:
    if (c == 2) {
      // L3C0 type
      fNrows        =  12;
      fLength       = 108.0;
      fWidth        = 105.5;
      fLengthOPad   =   8.0;
      fWidthOPad    =   0.775;
      fLengthIPad   =   9.0;
      fWidthIPad    =   0.725;
      fTiltingAngle =   2.0;
    }
    else {
      // L3C1 type
      fNrows        =  16;
      fLength       = 136.0;
      fWidth        = 105.5;
      fLengthOPad   =   7.5;
      fWidthOPad    =   0.775;
      fLengthIPad   =   8.5;
      fWidthIPad    =   0.725;
      fTiltingAngle =   2.0;
    }
    break;
  case 4:
    if (c == 2) {
      // L4C0 type
      fNrows        =  12;
      fLength       = 108.0;
      fWidth        = 109.9;
      fLengthOPad   =   8.0;
      fWidthOPad    =   0.845;
      fLengthIPad   =   9.0;
      fWidthIPad    =   0.755;
      fTiltingAngle =  -2.0;
    }
    else {
      // L4C1 type
      fNrows        =  16;
      fLength       = 143.0;
      fWidth        = 109.9;
      fLengthOPad   =   7.5;
      fWidthOPad    =   0.845;
      fLengthIPad   =   9.0;
      fWidthIPad    =   0.755;
      fTiltingAngle =  -2.0;
    }
    break;
  case 5:
    if (c == 2) {
      // L5C0 type
      fNrows        =  12;
      fLength       = 108.0;
      fWidth        = 114.4;
      fLengthOPad   =   8.0;
      fWidthOPad    =   0.965;
      fLengthIPad   =   9.0;
      fWidthIPad    =   0.785;
      fTiltingAngle =   2.0;
    }
    else {
      // L5C1 type
      fNrows        =  16;
      fLength       = 145.0;
      fWidth        = 114.4;
      fLengthOPad   =   8.5;
      fWidthOPad    =   0.965;
      fLengthIPad   =   9.0;
      fWidthIPad    =   0.785;
      fTiltingAngle =   2.0;
    }
    break;
  };

  //
  // The positions of the borders of the pads
  //
  // Row direction
  //
  if (fPadRow) delete [] fPadRow;
  fPadRow = new Double_t[fNrows];
  Double_t row = fGeo->GetChamberLength(p,0)
    	       + fGeo->GetChamberLength(p,1)
               + fGeo->GetChamberLength(p,2) / 2.
               - fLengthRim;
  for (Int_t ic = 0; ic < c; ic++) {
    row -= fGeo->GetChamberLength(p,ic);
  }
  for (Int_t ir = 0; ir < fNrows; ir++) {
    fPadRow[ir] = row;
    row -= fRowSpacing;
    if (ir == 1) {
      row -= fLengthOPad;
    }
    else {
      row -= fLengthIPad;
    }
  }
  //
  // Column direction
  //
  if (fPadCol) delete [] fPadCol;
  fPadCol = new Double_t[fNcols];
  Double_t col = fGeo->GetChamberWidth(p) / 2. 
               - fWidthRim;
  for (Int_t ic = 0; ic < fNcols; ic++) {
    fPadCol[ic] = col;
    col -= fColSpacing;
    if (ic == 1) {
      col -= fWidthOPad;
    }
    else {
      col -= fWidthIPad;
    }
  }

}

//_____________________________________________________________________________
AliTRDpadPlane::AliTRDpadPlane(const AliTRDpadPlane &p):TObject(p)
{
  //
  // AliTRDpadPlane copy constructor
  //

  ((AliTRDpadPlane &) p).Copy(*this);

}

//_____________________________________________________________________________
AliTRDpadPlane::~AliTRDpadPlane()
{
  //
  // AliTRDpadPlane destructor
  //

  if (fGeo) {
    delete fGeo;
    fGeo    = 0;
  }

  if (fPadRow) {
    delete [] fPadRow;
    fPadRow = 0;
  }

  if (fPadCol) {
    delete [] fPadCol;
    fPadCol = 0;
  }

}

//_____________________________________________________________________________
void AliTRDpadPlane::Copy(TObject &p) const
{
  //
  // Copy function
  //

  Int_t iBin = 0;

  ((AliTRDpadPlane &) p).fGeo          = 0;

  ((AliTRDpadPlane &) p).fPla          = fPla;
  ((AliTRDpadPlane &) p).fCha          = fCha;

  ((AliTRDpadPlane &) p).fLength       = fLength;
  ((AliTRDpadPlane &) p).fWidth        = fWidth;
  ((AliTRDpadPlane &) p).fLengthRim    = fLengthRim;
  ((AliTRDpadPlane &) p).fWidthRim     = fWidthRim;
  ((AliTRDpadPlane &) p).fLengthOPad   = fLengthOPad;
  ((AliTRDpadPlane &) p).fWidthOPad    = fWidthOPad;
  ((AliTRDpadPlane &) p).fLengthIPad   = fLengthIPad;
  ((AliTRDpadPlane &) p).fWidthIPad    = fWidthIPad;

  ((AliTRDpadPlane &) p).fRowSpacing   = fRowSpacing;
  ((AliTRDpadPlane &) p).fColSpacing   = fColSpacing;

  ((AliTRDpadPlane &) p).fNrows        = fNrows;
  ((AliTRDpadPlane &) p).fNcols        = fNcols;

  ((AliTRDpadPlane &) p).fTiltingAngle = fTiltingAngle;

  if (((AliTRDpadPlane &) p).fPadRow) delete [] ((AliTRDpadPlane &) p).fPadRow;
  ((AliTRDpadPlane &) p).fPadRow = new Double_t[fNrows];
  for (iBin = 0; iBin < fNrows; iBin++) {
    ((AliTRDpadPlane &) p).fPadRow[iBin] = fPadRow[iBin];
  }                                                                             

  if (((AliTRDpadPlane &) p).fPadCol) delete [] ((AliTRDpadPlane &) p).fPadCol;
  ((AliTRDpadPlane &) p).fPadCol = new Double_t[fNrows];
  for (iBin = 0; iBin < fNrows; iBin++) {
    ((AliTRDpadPlane &) p).fPadCol[iBin] = fPadCol[iBin];
  }                                                                             

  TObject::Copy(p);

}

//_____________________________________________________________________________
Int_t AliTRDpadPlane::GetPadRowNumber(const Double_t z)
{
  //
  // Finds the pad row number for a given global z-position
  //

  Int_t row    = 0;
  Int_t nabove = 0;
  Int_t nbelow = 0;
  Int_t middle = 0;

  if ((z > fPadRow[0]) ||
      (z < fPadRow[0] - fLength + 2.0*fLengthRim)) {
    row = -1;
  }
  else {
    nabove = fNrows+1;
    nbelow = 0;
    while (nabove - nbelow > 1) {
      middle = (nabove + nbelow) / 2;
      if (z == fPadRow[middle-1]) row    = middle;
      if (z  > fPadRow[middle-1]) nabove = middle;
      else                        nbelow = middle;
    }
    row = nbelow - 1;
  }

  return row;

}

//_____________________________________________________________________________
Int_t AliTRDpadPlane::GetPadColNumber(const Double_t rphi)
{
  //
  // Finds the pad column number for a given global rphi-position
  //

  Int_t col    = 0;
  Int_t nabove = 0;
  Int_t nbelow = 0;
  Int_t middle = 0;

  if ((rphi > fPadCol[0]) ||
      (rphi < fPadCol[0] - fWidth + 2.0*fWidthRim)) {
    col = -1;
  }
  else {
    nabove = fNcols+1;
    nbelow = 0;
    while (nabove - nbelow > 1) {
      middle = (nabove + nbelow) / 2;
      if (rphi == fPadCol[middle-1]) col    = middle;
      if (rphi  > fPadCol[middle-1]) nabove = middle;
      else                           nbelow = middle;
    }
    col = nbelow - 1;
  }

  return col;

}

//_____________________________________________________________________________
Double_t AliTRDpadPlane::GetPadRowOffset(const Int_t row, const Double_t z)
{
  //
  // Calculates the distance to the pad border in row direction
  //

  if ((row < 0) || (row >= fNrows)) {
    return -1.0;
  }
  else {
    return fPadRow[row] - z;
  }

}

//_____________________________________________________________________________
Double_t AliTRDpadPlane::GetPadColOffset(const Int_t col, const Double_t rphi)
{
  //
  // Calculates the distance to the pad border in column direction
  //

  if ((col < 0) || (col >= fNcols)) {
    return -1.0;
  }
  else {
    return fPadCol[col] - rphi;
  }

}
