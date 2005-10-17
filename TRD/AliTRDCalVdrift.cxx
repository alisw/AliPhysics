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
//  TRD calibration class for Vdrift                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalVdrift.h"
#include "AliTRDCalROCVdrift.h"

ClassImp(AliTRDCalVdrift)

//_____________________________________________________________________________
AliTRDCalVdrift::AliTRDCalVdrift():TNamed()
{
  //
  // AliTRDCalVdrift default constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fROCVdrift[idet] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalVdrift::AliTRDCalVdrift(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTRDCalVdrift constructor
  //

  for (Int_t isec = 0; isec < kNsect; isec++) {
    for (Int_t ipla = 0; ipla < kNplan; ipla++) {
      for (Int_t icha = 0; icha < kNcham; icha++) {
        Int_t idet = GetDet(ipla,icha,isec);
        fROCVdrift[idet] = new AliTRDCalROCVdrift(ipla,icha);
      }
    }
  }

}


//_____________________________________________________________________________
AliTRDCalVdrift::AliTRDCalVdrift(const AliTRDCalVdrift &c):TNamed(c)
{
  //
  // AliTRDCalVdrift copy constructor
  //

  ((AliTRDCalVdrift &) c).Copy(*this);

}

///_____________________________________________________________________________
AliTRDCalVdrift::~AliTRDCalVdrift()
{
  //
  // AliTRDCalVdrift destructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROCVdrift[idet]) {
      delete fROCVdrift[idet];
      fROCVdrift[idet] = 0;
    }
  }

}

//_____________________________________________________________________________
AliTRDCalVdrift &AliTRDCalVdrift::operator=(const AliTRDCalVdrift &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDCalVdrift &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDCalVdrift::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROCVdrift[idet]) {
      fROCVdrift[idet]->Copy(*((AliTRDCalVdrift &) c).fROCVdrift[idet]);
    }
  }

  TObject::Copy(c);

}

