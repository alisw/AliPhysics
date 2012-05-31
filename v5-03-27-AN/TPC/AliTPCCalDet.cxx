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
//  TPC calibration class for parameters which saved per detector            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCCalDet.h"

ClassImp(AliTPCCalDet)

//_____________________________________ ________________________________________
AliTPCCalDet::AliTPCCalDet():TNamed()
{
  //
  // AliTPCCalDet default constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fData[idet] = 0;
  }

}

//_____________________________________________________________________________
AliTPCCalDet::AliTPCCalDet(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTPCCalDet constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fData[idet] = 0;
  }

}


//_____________________________________________________________________________
AliTPCCalDet::AliTPCCalDet(const AliTPCCalDet &c):TNamed(c)
{
  //
  // AliTPCCalDet copy constructor
  //

  ((AliTPCCalDet &) c).Copy(*this);

}

///_____________________________________________________________________________
AliTPCCalDet::~AliTPCCalDet()
{
  //
  // AliTPCCalDet destructor
  //

}

//_____________________________________________________________________________
AliTPCCalDet &AliTPCCalDet::operator=(const AliTPCCalDet &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTPCCalDet &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTPCCalDet::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    ((AliTPCCalDet &) c).fData[idet] = fData[idet];
  }

  TObject::Copy(c);

}

