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

/* $Id: AliTRDUshortInfo.cxx 27946 2008-08-13 15:26:24Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Calibration base class for a single ROC                                  //
//  Contains one UShort_t value per pad                                      //
//  However, values are set and get as float, there are stored internally as //
//  (UShort_t) value * 10000                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDUshortInfo.h"

ClassImp(AliTRDUshortInfo)

//_____________________________________________________________________________
AliTRDUshortInfo::AliTRDUshortInfo()
  :TObject()
  ,fSize(0)
  ,fData(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDUshortInfo::AliTRDUshortInfo(Int_t n)
  :TObject()
  ,fSize(n)
  ,fData(0)
{
  //
  // Constructor that initializes a given size
  //
  
  fData = new UShort_t[fSize];
  for(Int_t k = 0; k < fSize; k++){
    fData[k] = 0;
  }

}

//_____________________________________________________________________________
AliTRDUshortInfo::AliTRDUshortInfo(const AliTRDUshortInfo &c)
  :TObject(c)
  ,fSize(c.fSize)
  ,fData(0)
{
  //
  // AliTRDUshortInfo copy constructor
  //

  Int_t iBin = 0;

  fData = new UShort_t[fSize];
  for (iBin = 0; iBin < fSize; iBin++) {
    fData[iBin] = ((AliTRDUshortInfo &) c).fData[iBin];
  }

}

//_____________________________________________________________________________
AliTRDUshortInfo::~AliTRDUshortInfo()
{
  //
  // AliTRDUshortInfo destructor
  //

  if (fData) {
    delete [] fData;
    fData = 0;
  }

}

//_____________________________________________________________________________
AliTRDUshortInfo &AliTRDUshortInfo::operator=(const AliTRDUshortInfo &c)
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
  fData = new UShort_t[fSize];
  for (Int_t iBin = 0; iBin < fSize; iBin++) {
    fData[iBin] = c.fData[iBin];
  }

  return *this;

}

//_____________________________________________________________________________
void AliTRDUshortInfo::Copy(TObject &c) const
{
  //
  // Copy function
  //
  
  Int_t iBin = 0;

  ((AliTRDUshortInfo &) c).fSize = fSize;

  if (((AliTRDUshortInfo &) c).fData) delete [] ((AliTRDUshortInfo &) c).fData;
  ((AliTRDUshortInfo &) c).fData = new UShort_t[fSize];
  for (iBin = 0; iBin < fSize; iBin++) {
    ((AliTRDUshortInfo &) c).fData[iBin] = fData[iBin];
  }
  
  TObject::Copy(c);

}

//_____________________________________________________________________________
void AliTRDUshortInfo::SetSize(Int_t n)
{
  //
  // Set the size
  //

  if(fData) delete [] fData;
  fData = new UShort_t[n];

  fSize = n;
  
}
