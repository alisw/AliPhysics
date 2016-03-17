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
#include <TBits.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <memory>

#include "AliEmcalTrackSelectionESD.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliVCuts.h"
#include "AliEmcalESDTrackCutsGenerator.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTrackSelectionESD)
/// \endcond

/**
 * Default constructor
 */
AliEmcalTrackSelectionESD::AliEmcalTrackSelectionESD():
		AliEmcalTrackSelection()
{
}

/**
 * Constructor with cuts
 */
AliEmcalTrackSelectionESD::AliEmcalTrackSelectionESD(AliVCuts* cuts):
		AliEmcalTrackSelection()
{
  this->AddTrackCuts(cuts);
}

/**
 * Constructor, initalising track cuts depending on the requested type of filtering
 *
 * \param type Track filtering type
 * \param period  Period string (e.g. LHC11h)
 */
AliEmcalTrackSelectionESD::AliEmcalTrackSelectionESD(ETrackFilterType_t type, const char* period):
  AliEmcalTrackSelection()
{
  GenerateTrackCuts(type, period);
}

/**
 * Automatically generates track cuts depending on the requested type of filtering
 *
 * \param type    Track filtering type
 * \param period  Period string (e.g. LHC11h)
 */
void AliEmcalTrackSelectionESD::GenerateTrackCuts(ETrackFilterType_t type, const char* period)
{
  if (fListOfCuts) fListOfCuts->Clear();
  fSelectionModeAny = kTRUE;

  switch (type) {
  case kHybridTracks:
    AliEmcalESDTrackCutsGenerator::AddHybridTrackCuts(this, period);
    break;

  case kTPCOnlyTracks:
    AliEmcalESDTrackCutsGenerator::AddTPCOnlyTrackCuts(this, period);
    break;

  default:
    break;
  }
}

/**
 * Check whether track is accepted. Iterates over all cuts assigned to the track selection.
 *
 * \param trk: Track to check
 * \return: true if selected, false otherwise
 */
bool AliEmcalTrackSelectionESD::IsTrackAccepted(AliVTrack* const trk) {
  if (!fListOfCuts) return kTRUE;
  AliESDtrack *esdt = dynamic_cast<AliESDtrack *>(trk);
  if (!esdt) {
    AliPicoTrack *picoTrack = dynamic_cast<AliPicoTrack *>(trk);
    if (picoTrack) {
      esdt = dynamic_cast<AliESDtrack*>(picoTrack->GetTrack());
    }
    else {
      AliError("Neither Pico nor ESD track");
      return kFALSE;
    }
  }

  fTrackBitmap.ResetAllBits();
  Int_t cutcounter = 0;
  for (TIter cutIter = TIter(fListOfCuts).Begin(); cutIter != TIter::End(); ++cutIter){
    if((static_cast<AliVCuts *>(*cutIter))->IsSelected(esdt)) fTrackBitmap.SetBitNumber(cutcounter);
    cutcounter++;
  }
  // In case of ANY at least one bit has to be set, while in case of ALL all bits have to be set
  if (fSelectionModeAny){
    return fTrackBitmap.CountBits() > 0 || cutcounter == 0;
  } else {
    return fTrackBitmap.CountBits() == cutcounter;
  }
}
