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
#include "AliEMCALTriggerRawPatch.h"

#include <iostream>

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerRawPatch)
/// \endcond

AliEMCALTriggerRawPatch::AliEMCALTriggerRawPatch():
  fBitMask(0),
  fCol0(-1),
  fRow0(-1),
  fSize(-1),
  fADC(0)
{
}

AliEMCALTriggerRawPatch::AliEMCALTriggerRawPatch(Int_t col0, Int_t row0, Int_t size, Double_t adc):
  fBitMask(0),
  fCol0(col0),
  fRow0(row0),
  fSize(size),
  fADC(adc)
{
}

void AliEMCALTriggerRawPatch::PrintStream(std::ostream &stream) const {
  stream << "Patch: Col[" << fCol0 << "], Row[" << fRow0 << "] with size " << fSize << " and ADC " << fADC;
}

bool AliEMCALTriggerRawPatch::operator ==(const AliEMCALTriggerRawPatch &other) const {
  return fRow0 == other.fRow0 && fCol0 == other.fCol0 && fBitMask == other.fBitMask;
}

bool AliEMCALTriggerRawPatch::operator <(const AliEMCALTriggerRawPatch &other) const {
  return fADC < other.fADC;
}

std::ostream &operator<<(std::ostream &stream, const AliEMCALTriggerRawPatch &patch){
  patch.PrintStream(stream);
  return stream;
}
