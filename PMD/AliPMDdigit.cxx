/***************************************************************************
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
//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store digits for ALICE-PMD                         //
//                                                     //
//-----------------------------------------------------//
#include "Riostream.h"
#include "Rtypes.h"
#include "AliPMDdigit.h"
#include <stdio.h>

ClassImp(AliPMDdigit)

AliPMDdigit::AliPMDdigit()
{
  // Default Constructor
  fTrNumber   = 0;
  fDet        = 0;
  fSMNumber   = 0;
  fCellNumber = 0;
  fADC        = 0.;
}

AliPMDdigit::AliPMDdigit(Int_t trnumber, Int_t det, Int_t smnumber, 
			 Int_t cellnumber, Float_t adc)
{
  // Constructor
  fTrNumber   = trnumber;
  fDet        = det;
  fSMNumber   = smnumber;
  fCellNumber = cellnumber;
  fADC        = adc;
}
AliPMDdigit::AliPMDdigit(const AliPMDdigit& pmddigit):TObject(pmddigit) {
  //Copy Constructor 
  if(&pmddigit == this) return;
  this->fTrNumber   = pmddigit.fTrNumber;
  this->fDet        = pmddigit.fDet;
  this->fSMNumber   = pmddigit.fSMNumber;
  this->fCellNumber = pmddigit.fCellNumber;
  this->fADC        = pmddigit.fADC;
  return;
}
AliPMDdigit & AliPMDdigit::operator=(const AliPMDdigit& pmddigit) {
  //Assignment operator 
  if(&pmddigit == this) return *this;
  this->fTrNumber   = pmddigit.fTrNumber;
  this->fDet        = pmddigit.fDet;
  this->fSMNumber   = pmddigit.fSMNumber;
  this->fCellNumber = pmddigit.fCellNumber;
  this->fADC        = pmddigit.fADC;
  return *this;
}
AliPMDdigit::~AliPMDdigit()
{
  // Default destructor
}
Int_t AliPMDdigit::GetTrackNumber() const
{
  return fTrNumber;
}
Int_t AliPMDdigit::GetDetector() const
{
  return fDet;
}
Int_t AliPMDdigit::GetSMNumber() const
{
  return fSMNumber;
}
Int_t AliPMDdigit::GetCellNumber() const
{
  return fCellNumber;
}
Float_t AliPMDdigit::GetADC() const
{
  return fADC;
}

