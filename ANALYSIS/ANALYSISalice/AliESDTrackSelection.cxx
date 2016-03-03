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
#include <AliESDTrackSelection.h>
#include <TBits.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <memory>

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliVCuts.h"

/// \cond CLASSIMP
ClassImp(AliESDTrackSelection)
/// \endcond

/**
 * Default constructor
 */
AliESDTrackSelection::AliESDTrackSelection():
		AliVTrackSelection()
{
}

/**
 * Constructor with cuts
 */
AliESDTrackSelection::AliESDTrackSelection(AliVCuts* cuts):
		AliVTrackSelection()
{
  this->AddTrackCuts(cuts);
}

/**
 * Check whether track is accepted. Iterates over all cuts assigned to the track selection.
 *
 * \param trk: Track to check
 * \return: true if selected, false otherwise
 */
bool AliESDTrackSelection::IsTrackAccepted(AliVTrack* const trk) {
  if (!fListOfCuts) return kTRUE;
  AliESDtrack *esdt = dynamic_cast<AliESDtrack *>(trk);
  if(!esdt){
    AliError("Failed getting ESD track");
    return kFALSE;
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
