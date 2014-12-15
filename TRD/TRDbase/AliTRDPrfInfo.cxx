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

/* $Id: AliTRDPrfInfo.cxx 27946 2008-08-13 15:26:24Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Calibration base class for a single ROC                                  //
//  Contains one UShort_t value per pad                                      //
//  However, values are set and get as float, there are stored internally as //
//  (UShort_t) value * 10000                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDPrfInfo.h"

ClassImp(AliTRDPrfInfo)

//_____________________________________________________________________________
AliTRDPrfInfo::AliTRDPrfInfo()
  :TObject()
  ,fSize(0)
  ,fData(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDPrfInfo::AliTRDPrfInfo(Int_t n)
  :TObject()
  ,fSize(n)
  ,fData(0)
{
  //
  // Constructor that initializes a given size
  //

  fData = new UChar_t[n];
  for(Int_t k = 0; k < fSize; k++){
    fData[k] = 0;
  }

}

//_____________________________________________________________________________
AliTRDPrfInfo::AliTRDPrfInfo(const AliTRDPrfInfo &c)
  :TObject(c)
  ,fSize(c.fSize)
  ,fData(0)
{
  //
  // AliTRDPrfInfo copy constructor
  //

  Int_t iBin = 0;

  fData = new UChar_t[fSize];
  for (iBin = 0; iBin < fSize; iBin++) {
    fData[iBin] = ((AliTRDPrfInfo &) c).fData[iBin];
  }

}

//_____________________________________________________________________________
AliTRDPrfInfo::~AliTRDPrfInfo()
{
  //
  // AliTRDPrfInfo destructor
  //

  if (fData) {
    delete [] fData;
    fData = 0;
  }

}

//_____________________________________________________________________________
AliTRDPrfInfo &AliTRDPrfInfo::operator=(const AliTRDPrfInfo &c)
{
  //
  // Assignment operator
  //

  if (this == &c) {
    return *this;
  }

  fSize = c.fSize;

  if (fData) {
    delete [] fData;
  }
  fData = new UChar_t[fSize];
  for (Int_t iBin = 0; iBin < fSize; iBin++) {
    fData[iBin] = ((AliTRDPrfInfo &) c).fData[iBin];
  }

  return *this;

}

//_____________________________________________________________________________
void AliTRDPrfInfo::Copy(TObject &c) const
{
  //
  // Copy function
  //

  Int_t iBin = 0;

  ((AliTRDPrfInfo &) c).fSize = fSize;

  if (((AliTRDPrfInfo &) c).fData) delete [] ((AliTRDPrfInfo &) c).fData;
  ((AliTRDPrfInfo &) c).fData = new UChar_t[fSize];
  for (iBin = 0; iBin < fSize; iBin++) {
    ((AliTRDPrfInfo &) c).fData[iBin] = fData[iBin];
  }
  
  TObject::Copy(c);

}

//_____________________________________________________________________________
void AliTRDPrfInfo::SetSize(Int_t n)
{
  //
  // Set the size
  //

  if (fData) delete [] fData;
  fData = new UChar_t[n];

  fSize = n;
  
}
