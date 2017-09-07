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

#include <TClonesArray.h>
#include <TBits.h>
#include <TObjArray.h>

#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <AliEmcalTrackSelectionAOD.h>
#include <AliESDtrack.h>
#include <AliPicoTrack.h>

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

  default:
    break;
  }
}

bool AliEmcalTrackSelectionAOD::IsTrackAccepted(AliVTrack * const trk)
{
  AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(trk);
  if (!aodt){
    AliPicoTrack *picotrack = dynamic_cast<AliPicoTrack*>(trk);
    if(picotrack) {
      aodt = dynamic_cast<AliAODTrack *>(picotrack->GetTrack());
    }
    else {
      AliError("Track neither AOD track nor pico track");
      return kFALSE;
    }
  }
  if(!aodt){
    AliError("Failed getting AOD track");
    return kFALSE;
  }

  fTrackBitmap.ResetAllBits();
  UInt_t cutcounter(0);
  if (fFilterBits) {
    if(aodt->TestFilterBit(fFilterBits)) fTrackBitmap.SetBitNumber(cutcounter);
    cutcounter++;
  }
  if (fFilterHybridTracks) {
    if (aodt->IsHybridGlobalConstrainedGlobal()) {
      // If the hybrid filter bits are not provided (fHybridFilterBits[0] == 0) all hybrid tracks will be selected in the same group
      if (fHybridFilterBits[0] < 0 || aodt->TestFilterBit(BIT(fHybridFilterBits[0]))) fTrackBitmap.SetBitNumber(cutcounter);
      if (aodt->TestFilterBit(BIT(fHybridFilterBits[1]))) fTrackBitmap.SetBitNumber(cutcounter+1);
    }
    cutcounter += 2;
  }
  if (fFilterTPCTracks) {
    if(aodt->IsHybridTPCConstrainedGlobal()) fTrackBitmap.SetBitNumber(cutcounter);
    cutcounter++;
  }
  if (fListOfCuts) {
    for (auto cutIter : *fListOfCuts){
      AliVCuts *trackCuts = static_cast<AliVCuts*>(static_cast<AliEmcalManagedObject *>(cutIter)->GetObject());
      if (trackCuts->IsA() == AliESDtrackCuts::Class()) {
        // If track cuts are AliESDtrackCuts, the track needs to be converted to an AliESDtrack before
        AliESDtrack copyTrack(aodt);
        if (trackCuts->IsSelected(&copyTrack)) fTrackBitmap.SetBitNumber(cutcounter);
      }
      else{
        if (trackCuts->IsSelected(aodt)) fTrackBitmap.SetBitNumber(cutcounter);
      }
      cutcounter++;
    }
  }

  if (fSelectionModeAny){
    // In case of ANY one of the cuts need to be fulfilled (equivalent to one but set)
    return fTrackBitmap.CountBits() > 0 || cutcounter == 0;
  }
  else {
    // In case of ALL all of the cuts need to be fulfilled (equivalent to all bits set)
    return fTrackBitmap.CountBits() == cutcounter;
  }
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
      (period.Length() == 6 && period.BeginsWith("lhc15")) // all Run-2 data, excluding MC productions
  ) {
    bits[0] = 8;
    bits[1] = 9;
  }

  else if (period == "lhc10f7a" || period == "lhc12a15e" || period.BeginsWith("lhc12a17") ||
      period == "lhc13b4" || period == "lhc13b4_fix" || period == "lhc13b4_plus" || period == "lhc14k1a" || period == "lhc14k1b" || period == "lhc13e4" ||
      period.BeginsWith("lhc14a1") || period.BeginsWith("lhc13b2_efix") ||
      period.BeginsWith("lhc15g6")) {
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
