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
//  TRD calibration class for the single pad status                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalPadStatus.h"
#include "AliTRDgeometry.h"
#include "AliTRDCalSingleChamberStatus.h"

ClassImp(AliTRDCalPadStatus)

//_____________________________________________________________________________
AliTRDCalPadStatus::AliTRDCalPadStatus()
  :TNamed()
{
  //
  // AliTRDCalPadStatus default constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fROC[idet] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalPadStatus::AliTRDCalPadStatus(const Text_t *name, const Text_t *title)
  :TNamed(name,title)
{
  //
  // AliTRDCalPadStatus constructor
  //

  for (Int_t isec = 0; isec < kNsect; isec++) {
    for (Int_t ipla = 0; ipla < kNplan; ipla++) {
      for (Int_t icha = 0; icha < kNcham; icha++) {
        Int_t idet = AliTRDgeometry::GetDetector(ipla,icha,isec);
        fROC[idet] = new AliTRDCalSingleChamberStatus(ipla,icha,144);
      }
    }
  }

}

//_____________________________________________________________________________
AliTRDCalPadStatus::AliTRDCalPadStatus(const AliTRDCalPadStatus &c)
  :TNamed(c)
{
  //
  // AliTRDCalPadStatus copy constructor
  //

  ((AliTRDCalPadStatus &) c).Copy(*this);

}

///_____________________________________________________________________________
AliTRDCalPadStatus::~AliTRDCalPadStatus()
{
  //
  // AliTRDCalPadStatus destructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]) {
      delete fROC[idet];
      fROC[idet] = 0;
    }
  }

}

//_____________________________________________________________________________
AliTRDCalPadStatus &AliTRDCalPadStatus::operator=(const AliTRDCalPadStatus &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDCalPadStatus &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDCalPadStatus::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]) {
      fROC[idet]->Copy(*((AliTRDCalPadStatus &) c).fROC[idet]);
    }
  }

  TObject::Copy(c);

}

//_____________________________________________________________________________
Bool_t AliTRDCalPadStatus::CheckStatus(Int_t d, Int_t col, Int_t row, Int_t bitMask) const
{
  //
  // Checks the pad status
  //

  AliTRDCalSingleChamberStatus *roc = GetCalROC(d);
  if (!roc) {
    return kFALSE;
  }
  else {
    return (roc->GetStatus(col, row) & bitMask) ? kTRUE : kFALSE;
  }

}

//_____________________________________________________________________________
AliTRDCalSingleChamberStatus* AliTRDCalPadStatus::GetCalROC(Int_t p, Int_t c, Int_t s) const
{ 
  //
  // Returns the readout chamber of this pad
  //

  return fROC[AliTRDgeometry::GetDetector(p,c,s)];   

}
