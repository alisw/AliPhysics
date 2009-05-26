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

AliPMDdigit::AliPMDdigit():
  fTrNumber(0),
  fTrPid(0),
  fDet(0),
  fSMNumber(0),
  fRow(0),
  fColumn(0),
  fADC(0.)
{
  // Default Constructor
}

AliPMDdigit::AliPMDdigit(Int_t trnumber, Int_t trpid, Int_t det,
			 Int_t smnumber, 
			 Int_t irow, Int_t icol, Float_t adc):
  fTrNumber(trnumber),
  fTrPid(trpid),
  fDet(det),
  fSMNumber(smnumber),
  fRow(irow),
  fColumn(icol),
  fADC(adc)
{
  // Constructor
}
AliPMDdigit::AliPMDdigit(AliPMDdigit *pmddigit):
  fTrNumber(0),
  fTrPid(0),
  fDet(0),
  fSMNumber(0),
  fRow(0),
  fColumn(0),
  fADC(0.)
{
  *this = *pmddigit;
}

AliPMDdigit::AliPMDdigit(const AliPMDdigit& pmddigit):
  TObject(pmddigit),
  fTrNumber(pmddigit.fTrNumber),
  fTrPid(pmddigit.fTrPid),
  fDet(pmddigit.fDet),
  fSMNumber(pmddigit.fSMNumber),
  fRow(pmddigit.fRow),
  fColumn(pmddigit.fColumn),
  fADC(pmddigit.fADC)
{
  //Copy Constructor 
}
AliPMDdigit & AliPMDdigit::operator=(const AliPMDdigit& pmddigit) {
  //Assignment operator 
  if(this != &pmddigit)
    {
      fTrNumber   = pmddigit.fTrNumber;
      fTrPid      = pmddigit.fTrPid;
      fDet        = pmddigit.fDet;
      fSMNumber   = pmddigit.fSMNumber;
      fRow        = pmddigit.fRow;
      fColumn     = pmddigit.fColumn;
      fADC        = pmddigit.fADC;
    }
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
Int_t AliPMDdigit::GetTrackPid() const
{
  return fTrPid;
}
Int_t AliPMDdigit::GetDetector() const
{
  return fDet;
}
Int_t AliPMDdigit::GetSMNumber() const
{
  return fSMNumber;
}
Int_t AliPMDdigit::GetRow() const
{
  return fRow;
}
Int_t AliPMDdigit::GetColumn() const
{
  return fColumn;
}
Float_t AliPMDdigit::GetADC() const
{
  return fADC;
}

