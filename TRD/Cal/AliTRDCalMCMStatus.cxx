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
//  TRD calibration class for MCM status                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalMCMStatus.h"
#include "AliTRDgeometry.h"
#include "AliTRDCalSingleChamberStatus.h"

ClassImp(AliTRDCalMCMStatus)

//_____________________________________________________________________________
AliTRDCalMCMStatus::AliTRDCalMCMStatus()
                   :TNamed()
{
  //
  // AliTRDCalMCMStatus default constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fROC[idet] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalMCMStatus::AliTRDCalMCMStatus(const Text_t *name, const Text_t *title)
                   :TNamed(name,title)
{
  //
  // AliTRDCalMCMStatus constructor
  //

  for (Int_t isec = 0; isec < kNsect; isec++) {
    for (Int_t ipla = 0; ipla < kNplan; ipla++) {
      for (Int_t icha = 0; icha < kNcham; icha++) {
        Int_t idet = AliTRDgeometry::GetDetector(ipla,icha,isec);
        fROC[idet] = new AliTRDCalSingleChamberStatus(ipla,icha,8);
      }
    }
  }

}

//_____________________________________________________________________________
AliTRDCalMCMStatus::AliTRDCalMCMStatus(const AliTRDCalMCMStatus &c)
                   :TNamed(c)
{
  //
  // AliTRDCalMCMStatus copy constructor
  //

  ((AliTRDCalMCMStatus &) c).Copy(*this);

}

///_____________________________________________________________________________
AliTRDCalMCMStatus::~AliTRDCalMCMStatus()
{
  //
  // AliTRDCalMCMStatus destructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]) {
      delete fROC[idet];
      fROC[idet] = 0;
    }
  }

}

//_____________________________________________________________________________
AliTRDCalMCMStatus &AliTRDCalMCMStatus::operator=(const AliTRDCalMCMStatus &c)
{
  //
  // Assignment operator
  //

  if (this != &c) {
    ((AliTRDCalMCMStatus &) c).Copy(*this);
  }

  return *this;

}

//_____________________________________________________________________________
void AliTRDCalMCMStatus::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]) {
      fROC[idet]->Copy(*((AliTRDCalMCMStatus &) c).fROC[idet]);
    }
  }

  TObject::Copy(c);

}

//_____________________________________________________________________________
Bool_t AliTRDCalMCMStatus::CheckStatus(Int_t d, Int_t col, Int_t row, Int_t bitMask) const
{
  //
  // Checks the MCM status byte
  //

  // To translate pad column number into MCM number
  Int_t mcm = ((Int_t) col / 18);

  AliTRDCalSingleChamberStatus* roc = GetCalROC(d);

  if (!roc) {
    return kFALSE;
  }
  else {
    return (roc->GetStatus(mcm,row) & bitMask) ? kTRUE : kFALSE;
  }

}

//_____________________________________________________________________________
Char_t AliTRDCalMCMStatus::GetStatus(Int_t d, Int_t col, Int_t row) const
{
  //
  // Gets the MCM status byte
  //

  // To translate pad column number into MCM number
  Int_t mcm = ((Int_t) col / 18);

  AliTRDCalSingleChamberStatus* roc = GetCalROC(d);

  if (!roc) {
    return kFALSE;
  }
  else {
    return roc->GetStatus(mcm,row);
  }

}

//_____________________________________________________________________________
AliTRDCalSingleChamberStatus *AliTRDCalMCMStatus::GetCalROC(Int_t p, Int_t c, Int_t s) const
{
  //
  // Returns the readout chamber of this MCM
  //
 
  return GetCalROC(AliTRDgeometry::GetDetector(p,c,s)); 

}
