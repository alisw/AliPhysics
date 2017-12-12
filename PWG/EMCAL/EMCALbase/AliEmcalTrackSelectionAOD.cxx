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
#include <vector>

#include <TClonesArray.h>
#include <TBits.h>
#include <TObjArray.h>

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliEmcalAODFilterBitCuts.h"
#include "AliEmcalAODHybridTrackCuts.h"
#include "AliEmcalCutBase.h"
#include "AliEmcalTrackSelResultCombined.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliPicoTrack.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTrackSelectionAOD)
/// \endcond

AliEmcalTrackSelectionAOD::AliEmcalTrackSelectionAOD() :
	AliEmcalTrackSelection(),
	fFilterBits(0),
	fFilterHybridTracks(kFALSE),
	fFilterTPCTracks(kFALSE)
{
  fHybridFilterBits[0] = -1;
  fHybridFilterBits[1] = -1;
}

AliEmcalTrackSelectionAOD::AliEmcalTrackSelectionAOD(AliVCuts* cuts, UInt_t filterbits):
	AliEmcalTrackSelection(),
	fFilterBits(filterbits),
  fFilterHybridTracks(kFALSE),
  fFilterTPCTracks(kFALSE)
{
  fHybridFilterBits[0] = -1;
  fHybridFilterBits[1] = -1;
  if(cuts) AddTrackCuts(cuts);
}

AliEmcalTrackSelectionAOD::AliEmcalTrackSelectionAOD(ETrackFilterType_t type, const char* period):
  AliEmcalTrackSelection(),
  fFilterBits(0),
  fFilterHybridTracks(kFALSE),
  fFilterTPCTracks(kFALSE)
{
  fHybridFilterBits[0] = -1;
  fHybridFilterBits[1] = -1;
  GenerateTrackCuts(type, period);
}

void AliEmcalTrackSelectionAOD::GenerateTrackCuts(ETrackFilterType_t type, const char* period)
{
  switch (type) {
  case kHybridTracks:
    if (fListOfCuts) fListOfCuts->Clear();
    fFilterBits = 0;
    fFilterHybridTracks = kTRUE;
    fFilterTPCTracks = kFALSE;
    GetHybridFilterBits(fHybridFilterBits, period);
    fSelectionModeAny = kTRUE;
    break;

  case kTPCOnlyTracks:
    if (fListOfCuts) fListOfCuts->Clear();
    fFilterBits = 0;
    fFilterHybridTracks = kFALSE;
    fHybridFilterBits[0] = -1;
    fHybridFilterBits[1] = -1;
    fFilterTPCTracks = kTRUE;
    break;

  case AliEmcalTrackSelection::kITSPureTracks:
    if (fListOfCuts) fListOfCuts->Clear();
    AddTrackCuts(AliESDtrackCuts::GetStandardITSSATrackCuts2010());
    fFilterHybridTracks = kFALSE;
    fFilterTPCTracks = kFALSE;
    fHybridFilterBits[0] = -1;
    fHybridFilterBits[1] = -1;
    break;

  case kHybridTracks2010wNoRefit:
    {
      auto trackcuts = new PWG::EMCAL::AliEmcalAODHybridTrackCuts("hybridcuts2010_wNoRefit");
      AddTrackCuts(trackcuts);
      break;
    }

  case kHybridTracks2010woNoRefit:
    {
      auto trackcuts = new PWG::EMCAL::AliEmcalAODHybridTrackCuts("hybridcuts2010_woNoRefit");
      trackcuts->SetSelectNonITSrefitTracks(kFALSE);
      AddTrackCuts(trackcuts);
      break;
    }

  case kHybridTracks2011wNoRefit:
    {
      auto trackcuts = new PWG::EMCAL::AliEmcalAODHybridTrackCuts("hybridcuts2011_wNoRefit");
      AddTrackCuts(trackcuts);
      break;
    }

  case kHybridTracks2011woNoRefit:
    {
      auto trackcuts = new PWG::EMCAL::AliEmcalAODHybridTrackCuts("hybridcuts2011_woNoRefit");
      trackcuts->SetSelectNonITSrefitTracks(kFALSE);
      AddTrackCuts(trackcuts);
      break;
    }

  default:
    break;
  }
}

PWG::EMCAL::AliEmcalTrackSelResultPtr AliEmcalTrackSelectionAOD::IsTrackAccepted(AliVTrack * const trk)
{
  AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(trk);
  if (!aodt){
    AliPicoTrack *picotrack = dynamic_cast<AliPicoTrack*>(trk);
    if(picotrack) {
      aodt = dynamic_cast<AliAODTrack *>(picotrack->GetTrack());
    }
    else {
      AliError("Track neither AOD track nor pico track");
      return PWG::EMCAL::AliEmcalTrackSelResultPtr(nullptr, kFALSE);
    }
  }
  if(!aodt){
    AliError("Failed getting AOD track");
    return PWG::EMCAL::AliEmcalTrackSelResultPtr(nullptr, kFALSE);
  }

  TBits trackbitmap(64);
  trackbitmap.ResetAllBits();
  UInt_t cutcounter(0);
  std::vector<PWG::EMCAL::AliEmcalTrackSelResultPtr> selectionStatus;
  if (fFilterBits) {
    if(aodt->TestFilterBit(fFilterBits)) {
      trackbitmap.SetBitNumber(cutcounter);
      selectionStatus.emplace_back(PWG::EMCAL::AliEmcalTrackSelResultPtr(aodt, kTRUE));
    } else {
      selectionStatus.emplace_back(PWG::EMCAL::AliEmcalTrackSelResultPtr(aodt, kFALSE));
    }
    cutcounter++;
  }
  if (fFilterHybridTracks) {
    if (aodt->IsHybridGlobalConstrainedGlobal()) {
      // If the hybrid filter bits are not provided (fHybridFilterBits[0] == 0) all hybrid tracks will be selected in the same group
      if (fHybridFilterBits[0] < 0 || aodt->TestFilterBit(BIT(fHybridFilterBits[0]))) {
        selectionStatus.emplace_back(PWG::EMCAL::AliEmcalTrackSelResultPtr(aodt, kTRUE));
        trackbitmap.SetBitNumber(cutcounter);
      } else {
        selectionStatus.emplace_back(PWG::EMCAL::AliEmcalTrackSelResultPtr(aodt, kFALSE));
      } 
      if (aodt->TestFilterBit(BIT(fHybridFilterBits[1]))) trackbitmap.SetBitNumber(cutcounter+1);
    }
    cutcounter += 2;
  }
  if (fFilterTPCTracks) {
    if(aodt->IsHybridTPCConstrainedGlobal()) {
      selectionStatus.emplace_back(PWG::EMCAL::AliEmcalTrackSelResultPtr(aodt, kTRUE));
      trackbitmap.SetBitNumber(cutcounter);
    } else {
      selectionStatus.emplace_back(PWG::EMCAL::AliEmcalTrackSelResultPtr(aodt, kTRUE));
    }
    cutcounter++;
  }
  if (fListOfCuts) {
    for (auto cutIter : *fListOfCuts){
      PWG::EMCAL::AliEmcalCutBase *trackCuts = static_cast<PWG::EMCAL::AliEmcalCutBase*>(static_cast<AliEmcalManagedObject *>(cutIter)->GetObject());
      PWG::EMCAL::AliEmcalTrackSelResultPtr cutresults = trackCuts->IsSelected(aodt);
      if (cutresults) trackbitmap.SetBitNumber(cutcounter);
      selectionStatus.emplace_back(cutresults);
      cutcounter++;
    }
  }

  PWG::EMCAL::AliEmcalTrackSelResultPtr result(aodt, kFALSE, new PWG::EMCAL::AliEmcalTrackSelResultCombined(selectionStatus));
  if (fSelectionModeAny){
    // In case of ANY one of the cuts need to be fulfilled (equivalent to one but set)
    result.SetSelectionResult(trackbitmap.CountBits() > 0 || cutcounter == 0);
  }
  else {
    // In case of ALL all of the cuts need to be fulfilled (equivalent to all bits set)
    result.SetSelectionResult(trackbitmap.CountBits() == cutcounter);
  }
  return result;
}

Bool_t AliEmcalTrackSelectionAOD::GetHybridFilterBits(Char_t bits[], TString period)
{
  period.ToLower();
  if (period == "lhc10b" || period == "lhc10c" || period == "lhc10d" ||
      period == "lhc10e" || period == "lhc10h" ||
      period == "lhc11h" || period == "lhc12a" || period == "lhc12b" ||
      period == "lhc12c" || period == "lhc12d" || period == "lhc12e" ||
      period == "lhc12f" || period == "lhc12g" || period == "lhc12h" ||
      period == "lhc12i" || period == "lhc13b" || period == "lhc13c" ||
      period == "lhc13d" || period == "lhc13e" || period == "lhc13f" ||
      period == "lhc13g" ||
      (period.Length() == 6 && (period.BeginsWith("lhc15") || period.BeginsWith("lhc16") || period.BeginsWith("lhc17"))) // all Run-2 data, excluding MC productions
  ) {
    bits[0] = 8;
    bits[1] = 9;
  }

  else if (period == "lhc10f7a" || period == "lhc12a15e" || period.BeginsWith("lhc12a17") ||
      period == "lhc13b4" || period == "lhc13b4_fix" || period == "lhc13b4_plus" || period == "lhc14k1a" || period == "lhc14k1b" || period == "lhc13e4" ||
      period.BeginsWith("lhc14a1") || period.BeginsWith("lhc13b2_efix") ||
      period.BeginsWith("lhc15g6") || period.BeginsWith("lhc16e1") || period.BeginsWith("lhc17f8")) {
    bits[0] = 8;
    bits[1] = 9;
  }

  else if (period == "lhc11a" || period == "lhc10hold" || period == "lhc11c" || period == "lhc11d") {
    bits[0] = 8;
    bits[1] = 4;
  }

  else if (period.Contains("lhc12a15a") || period == "lhc12a15f" ||
      period == "lhc12a15g" || period.BeginsWith("lhc11a1")) {
    bits[0] = 8;
    bits[1] = 4;
  }

  else {
    ::Error("AliEmcalTrackSelectionAOD::GetHybridFilterBits", "Could not find period %s! Hybrid tracks will be selected, but will not be able to distinguish between global and constrained.", period.Data());
    bits[0] = -1;
    bits[1] = -1;
    return kFALSE;
  }

  return kTRUE;
}
