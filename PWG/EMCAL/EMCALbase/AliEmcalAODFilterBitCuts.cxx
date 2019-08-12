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
#include <bitset>
#include <iostream>
#include "AliAODTrack.h"
#include "AliEmcalAODFilterBitCuts.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(PWG::EMCAL::AliEmcalAODFilterBitCuts)
/// \endcond

namespace PWG{

namespace EMCAL{

AliEmcalAODFilterBitCuts::AliEmcalAODFilterBitCuts():
  AliVCuts(),
  fAODfilterBits(0),
  fAODstatusBits(0),
  fSelectionMode(kSelAny)
{
}

AliEmcalAODFilterBitCuts::AliEmcalAODFilterBitCuts(const char *name, const char *title):
    AliVCuts(name, title),
    fAODfilterBits(0),
    fAODstatusBits(0),
    fSelectionMode(kSelAny)
{
}

Bool_t AliEmcalAODFilterBitCuts::IsSelected(TObject *o){
  AliDebugStream(1) << "Filter bit cut: selecting " << std::bitset<sizeof(decltype(fAODfilterBits))*8>(fAODfilterBits) << std::endl;
  AliAODTrack *testtrack = dynamic_cast<AliAODTrack *>(o);
  if(!testtrack) return false;
  Bool_t result(true);
  if(fAODfilterBits){
    if(!IsFilterBitsSelected(testtrack)) result = false;
  }
  if(fAODstatusBits){
    if(!IsStatusBitsSelected(testtrack)) result = false;
  }
  return result;
}

Bool_t AliEmcalAODFilterBitCuts::IsFilterBitsSelected(const AliAODTrack *const trk) const {
  Bool_t result(false);
  switch(fSelectionMode){
  case kSelAny: result = ((trk->GetFilterMap() & fAODfilterBits) > 0); break;
  case kSelAll: result = ((trk->GetFilterMap() & fAODfilterBits) == fAODfilterBits); break;
  };
  return result;
}

Bool_t AliEmcalAODFilterBitCuts::IsStatusBitsSelected(const AliAODTrack *const trk) const {
  ULong_t trackstatusbits = trk->TestBits(fAODstatusBits);
  Bool_t result(false);
  switch(fSelectionMode) {
  case kSelAny: result = (trackstatusbits > fAODstatusBits); break;
  case kSelAll: result = (trackstatusbits == fAODstatusBits); break;
  }
  return result;
}

}

}
