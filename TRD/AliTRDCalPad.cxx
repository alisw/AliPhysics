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
//  TRD calibration class for parameters which saved per pad                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalPad.h"
#include "AliTRDCalROC.h"
#include "AliTRDCalDet.h"

ClassImp(AliTRDCalPad)

//_____________________________________________________________________________
AliTRDCalPad::AliTRDCalPad():TNamed()
{
  //
  // AliTRDCalPad default constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fROC[idet] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalPad::AliTRDCalPad(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTRDCalPad constructor
  //

  for (Int_t isec = 0; isec < kNsect; isec++) {
    for (Int_t ipla = 0; ipla < kNplan; ipla++) {
      for (Int_t icha = 0; icha < kNcham; icha++) {
        Int_t idet = GetDet(ipla,icha,isec);
        fROC[idet] = new AliTRDCalROC(ipla,icha);
      }
    }
  }

}


//_____________________________________________________________________________
AliTRDCalPad::AliTRDCalPad(const AliTRDCalPad &c):TNamed(c)
{
  //
  // AliTRDCalPad copy constructor
  //

  ((AliTRDCalPad &) c).Copy(*this);

}

///_____________________________________________________________________________
AliTRDCalPad::~AliTRDCalPad()
{
  //
  // AliTRDCalPad destructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]) {
      delete fROC[idet];
      fROC[idet] = 0;
    }
  }

}

//_____________________________________________________________________________
AliTRDCalPad &AliTRDCalPad::operator=(const AliTRDCalPad &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDCalPad &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDCalPad::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]) {
      fROC[idet]->Copy(*((AliTRDCalPad &) c).fROC[idet]);
    }
  }

  TObject::Copy(c);
}

//_____________________________________________________________________________
void AliTRDCalPad::ScaleROCs(AliTRDCalDet* values)
{
  // 
  // Scales ROCs of this class with the values from the class <values>
  // Is used if an AliTRDCalPad object defines local variations of a parameter
  // defined per detector using a AliTRDCalDet class
  //
  
  if (!values)
    return;
  
  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]) { 
      fROC[idet]->Scale(values->GetValue(idet));
    }
  }
}

