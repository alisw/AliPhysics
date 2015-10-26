/**
 * @file AliEmcalTriggerRawPatchAP.cxx
 * @since Oct 23, 2015
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliEmcalTriggerRawPatchAP.h"
#include <iostream>

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerRawPatchAP)
/// \endcond

AliEmcalTriggerRawPatchAP::AliEmcalTriggerRawPatchAP():
  fBitMask(0),
  fCol0(-1),
  fRow0(-1),
  fSize(-1),
  fADC(0)
{
}


AliEmcalTriggerRawPatchAP::AliEmcalTriggerRawPatchAP(Int_t col0, Int_t row0, Int_t size, Double_t adc):
  fBitMask(0),
  fCol0(col0),
  fRow0(row0),
  fSize(size),
  fADC(adc)
{
}

void AliEmcalTriggerRawPatchAP::PrintStream(std::ostream &stream) const {
  stream << "Patch: Col[" << fCol0 << "], Row[" << fRow0 << "] with size " << fSize << " and ADC " << fADC;
}

bool AliEmcalTriggerRawPatchAP::operator ==(const AliEmcalTriggerRawPatchAP &other) const {
  return fRow0 == other.fRow0 && fCol0 == other.fCol0 && fBitMask == other.fBitMask;
}

bool AliEmcalTriggerRawPatchAP::operator <(const AliEmcalTriggerRawPatchAP &other) const {
  return fADC < other.fADC;
}

std::ostream &operator<<(std::ostream &stream, const AliEmcalTriggerRawPatchAP &patch){
  patch.PrintStream(stream);
  return stream;
}
