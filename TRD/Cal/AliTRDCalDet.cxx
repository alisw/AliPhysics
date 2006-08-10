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
//  TRD calibration class for parameters which saved per detector            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalDet.h"

ClassImp(AliTRDCalDet)

//_____________________________________________________________________________
AliTRDCalDet::AliTRDCalDet():TNamed()
{
  //
  // AliTRDCalDet default constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fData[idet] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalDet::AliTRDCalDet(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTRDCalDet constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fData[idet] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalDet::AliTRDCalDet(const AliTRDCalDet &c):TNamed(c)
{
  //
  // AliTRDCalDet copy constructor
  //

  ((AliTRDCalDet &) c).Copy(*this);

}

///_____________________________________________________________________________
AliTRDCalDet::~AliTRDCalDet()
{
  //
  // AliTRDCalDet destructor
  //

}

//_____________________________________________________________________________
AliTRDCalDet &AliTRDCalDet::operator=(const AliTRDCalDet &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDCalDet &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDCalDet::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    ((AliTRDCalDet &) c).fData[idet] = fData[idet];
  }

  TObject::Copy(c);

}

