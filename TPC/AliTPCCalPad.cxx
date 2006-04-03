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
//  TPC calibration class for parameters which saved per pad                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"

ClassImp(AliTPCCalPad)

//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad():TNamed()
{
  //
  // AliTPCCalPad default constructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    fROC[isec] = 0;
  }

}

//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTPCCalPad constructor
  //
  for (Int_t isec = 0; isec < kNsec; isec++) {
    fROC[isec] = new AliTPCCalROC(isec);
  }
}


//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad(const AliTPCCalPad &c):TNamed(c)
{
  //
  // AliTPCCalPad copy constructor
  //

  ((AliTPCCalPad &) c).Copy(*this);

}

///_____________________________________________________________________________
AliTPCCalPad::~AliTPCCalPad()
{
  //
  // AliTPCCalPad destructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    if (fROC[isec]) {
      delete fROC[isec];
      fROC[isec] = 0;
    }
  }

}

//_____________________________________________________________________________
AliTPCCalPad &AliTPCCalPad::operator=(const AliTPCCalPad &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTPCCalPad &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTPCCalPad::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    if (fROC[isec]) {
      fROC[isec]->Copy(*((AliTPCCalPad &) c).fROC[isec]);
    }
  }
  TObject::Copy(c);
}
