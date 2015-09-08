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
 * Select tracks from a TClonesArray of input tracks
 *
 * \param tracks TClonesArray of tracks (must not be null)
 * \return: TObjArray of selected tracks
 */
TObjArray* AliEmcalTrackSelectionESD::GetAcceptedTracks(const TClonesArray* const tracks) {
	if(!fListOfTracks) fListOfTracks = new TObjArray;
	else fListOfTracks->Clear();
	for(TIter trackIter = TIter(tracks).Begin(); trackIter != TIter::End(); ++trackIter){
	  if(IsTrackAccepted(static_cast<AliVTrack *>(*trackIter))) fListOfTracks->Add(*trackIter);
	}
	return fListOfTracks;
}

/**
 * Select tracks from a virtual event. Delegates selection process to function IsTrackAccepted
 *
 * \param event AliESDEvent, via interface of virtual event (must not be null)
 * \return TObjArray of selected tracks
 */
TObjArray* AliEmcalTrackSelectionESD::GetAcceptedTracks(const AliVEvent* const event) {
	if(!fListOfTracks) fListOfTracks = new TObjArray;
	else fListOfTracks->Clear();
	const AliESDEvent *esd = dynamic_cast<const AliESDEvent *>(event);
	if(!esd){
		AliError("Event not of type AliESDEvent");
		return fListOfTracks;
	}
	for(int itrk = 0; itrk < esd->GetNumberOfTracks(); itrk++){
	  AliESDtrack *trk = static_cast<AliESDtrack *>(esd->GetTrack(itrk));
	  if(IsTrackAccepted(trk)) fListOfTracks->AddLast(trk);
	}
	return fListOfTracks;
}

/**
 * Check whether track is accepted. Itterates over all cuts assinged to the track selection.
 *
 * \param trk: Track to check
 * \return: true if selected, false otherwise
 */
bool AliEmcalTrackSelectionESD::IsTrackAccepted(AliVTrack* const trk) {
  if(!fListOfCuts) return kTRUE;
  AliESDtrack *esdt = dynamic_cast<AliESDtrack *>(trk);
  if(!esdt){
    AliPicoTrack *picoTrack = dynamic_cast<AliPicoTrack *>(trk);
    if(picoTrack)
      esdt = dynamic_cast<AliESDtrack*>(picoTrack->GetTrack());
    else{
      AliError("Neither Pico nor ESD track");
      return kFALSE;
    }
  }

  TBits selectedMap(64);        // Counting track cuts among which track was SELECTED
  Int_t cutcounter = 0;
  for(TIter cutIter = TIter(fListOfCuts).Begin(); cutIter != TIter::End(); ++cutIter){
    if((static_cast<AliVCuts *>(*cutIter))->IsSelected(esdt)) selectedMap.SetBitNumber(cutcounter++);
  }
  // In case of ANY at least one bit has to be set, while in case of ALL all bits have to be set
  if(fSelectionModeAny){
    return selectedMap.CountBits(0) > 0;
  } else {
    return selectedMap.CountBits(0) == fListOfCuts->GetEntries();
  }
}
