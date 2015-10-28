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
#include "AliEmcalTriggerRawPatch.h"

#include <iostream>

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerRawPatch)
/// \endcond

/**
 * Dummy constructor
 */
AliEmcalTriggerRawPatch::AliEmcalTriggerRawPatch():
  fBitMask(0),
  fCol0(-1),
  fRow0(-1),
  fSize(-1),
  fADC(0)
{
}

/**
 * Main constructor
 * @param col0 Starting column
 * @param row0 Starting row
 * @param size Patch size
 * @param adc ADC value
 */
AliEmcalTriggerRawPatch::AliEmcalTriggerRawPatch(Int_t col0, Int_t row0, Int_t size, Double_t adc):
  fBitMask(0),
  fCol0(col0),
  fRow0(row0),
  fSize(size),
  fADC(adc)
{
}

/**
 * Print trigger patch information to a stream
 * @param stream Output stream
 */
void AliEmcalTriggerRawPatch::PrintStream(std::ostream &stream) const {
  stream << "Patch: Col[" << fCol0 << "], Row[" << fRow0 << "] with size " << fSize << " and ADC " << fADC;
}

/**
 * Comparison operator for equalness: Patches are equal if they have the same position and
 * the same trigger bit mask.
 * @param other Patch to compare to
 * @return True if the patches share the same position and trigger bit mask, false otherwise
 */
bool AliEmcalTriggerRawPatch::operator ==(const AliEmcalTriggerRawPatch &other) const {
  return fRow0 == other.fRow0 && fCol0 == other.fCol0 && fBitMask == other.fBitMask;
}

/**
 * Comparison operator for smaller. As this is used in sorting algorithms, the comparison
 * is made based on the patch ADC.
 * @param other Patch to compate to
 * @return True if the patch ADC of this patch is smaller, false otherwise
 */
bool AliEmcalTriggerRawPatch::operator <(const AliEmcalTriggerRawPatch &other) const {
  return fADC < other.fADC;
}

/**
 * output stream operator
 */
std::ostream &operator<<(std::ostream &stream, const AliEmcalTriggerRawPatch &patch){
  patch.PrintStream(stream);
  return stream;
}
