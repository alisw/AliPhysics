/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include <AliEmcalAODFilterBitCuts.h>
#include "AliAODTrack.h"

/// \cond CLASSIMP
ClassImp(AliEmcalAODFilterBitCuts)
/// \endcond

AliEmcalAODFilterBitCuts::AliEmcalAODFilterBitCuts():
  AliVCuts(),
  fAODfilterBits(0),
  fSelectionMode(kSelAny)
{
}

Bool_t AliEmcalAODFilterBitCuts::IsSelected(TObject *o){
  AliAODTrack *testtrack = dynamic_cast<AliAODTrack *>(o);
  if(!testtrack) return false;
  Bool_t result = false;
  switch(fSelectionMode){
  case kSelAny: result = ((testtrack->GetFilterMap() & fAODfilterBits) > 0); break;
  case kSelAll: result = ((testtrack->GetFilterMap() & fAODfilterBits) == fAODfilterBits);
  }
  return result;
}
