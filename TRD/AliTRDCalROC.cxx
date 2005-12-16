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
//  Calibration base class for a single ROC                                  //
//  Contains one float value per pad                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalROC.h"

ClassImp(AliTRDCalROC)

//_____________________________________________________________________________
AliTRDCalROC::AliTRDCalROC():TObject()
{
  //
  // Default constructor
  //

  fPla          = 0;
  fCha          = 0;

  fNrows        = 0;
  fNcols        = 0;

  fNchannels    = 0;
  fData         = 0;
}

//_____________________________________________________________________________
AliTRDCalROC::AliTRDCalROC(Int_t p, Int_t c):TObject()
{
  //
  // Constructor that initializes a given pad plane type
  //

  fPla = p;
  fCha = c;

  fNcols      = 144;

  //
  // The pad plane parameter
  //
  switch (p) {
  case 0:
    if (c == 2) {
      // L0C0 type
      fNrows        =  12;
    }
    else {
      // L0C1 type
      fNrows        =  16;
    }
    break;
  case 1:
    if (c == 2) {
      // L1C0 type
      fNrows        =  12;
    }
    else {
      // L1C1 type
      fNrows        =  16;
    }
    break;
  case 2:
    if (c == 2) {
      // L2C0 type
      fNrows        =  12;
    }
    else {
      // L2C1 type
      fNrows        =  16;
    }
    break;
  case 3:
    if (c == 2) {
      // L3C0 type
      fNrows        =  12;
    }
    else {
      // L3C1 type
      fNrows        =  16;
    }
    break;
  case 4:
    if (c == 2) {
      // L4C0 type
      fNrows        =  12;
    }
    else {
      // L4C1 type
      fNrows        =  16;
    }
    break;
  case 5:
    if (c == 2) {
      // L5C0 type
      fNrows        =  12;
    }
    else {
      // L5C1 type
      fNrows        =  16;
    }
    break;
  };

  fNchannels = fNrows * fNcols;
  if (fNchannels != 0)
    fData = new Float_t[fNchannels];
}

//_____________________________________________________________________________
AliTRDCalROC::AliTRDCalROC(const AliTRDCalROC &c):TObject(c)
{
  //
  // AliTRDCalROC copy constructor
  //

  ((AliTRDCalROC &) c).Copy(*this);

}

//_____________________________________________________________________________
AliTRDCalROC::~AliTRDCalROC()
{
  //
  // AliTRDCalROC destructor
  //

  if (fData) {
    delete [] fData;
    fData = 0;
  }
}

//_____________________________________________________________________________
AliTRDCalROC &AliTRDCalROC::operator=(const AliTRDCalROC &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDCalROC &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDCalROC::Copy(TObject &c) const
{
  //
  // Copy function
  //

  ((AliTRDCalROC &) c).fPla          = fPla;
  ((AliTRDCalROC &) c).fCha          = fCha;

  ((AliTRDCalROC &) c).fNrows        = fNrows;
  ((AliTRDCalROC &) c).fNcols        = fNcols;

  Int_t iBin = 0;

  ((AliTRDCalROC &) c).fNchannels = fNchannels;

  if (((AliTRDCalROC &) c).fData) delete [] ((AliTRDCalROC &) c).fData;
  ((AliTRDCalROC &) c).fData = new Float_t[fNchannels];
  for (iBin = 0; iBin < fNchannels; iBin++) {
    ((AliTRDCalROC &) c).fData[iBin] = fData[iBin];
  }

  TObject::Copy(c);

}

//_____________________________________________________________________________
void AliTRDCalROC::Scale(Float_t value)
{
  //
  // Scales all values of this ROC with the provided parameter. Is used if ROC defines
  // local variations of a global (or per detector defined) parameter
  //

  for (Int_t iBin = 0; iBin < fNchannels; iBin++) {
    fData[iBin] *= value;
  }
}
