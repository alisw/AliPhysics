/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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
/*
 * Implementation of the track selection for the analysis on ESDs using
 * AliESDtrackCuts as underlying structure
 *
 * Author:
 * 		Markus Fasel
 */
#include <TClonesArray.h>
#include <memory>

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"

#include <AliEMCalPtTaskTrackSelectionESD.h>

ClassImp(EMCalTriggerPtAnalysis::AliEMCalPtTaskTrackSelectionESD)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalPtTaskTrackSelectionESD::AliEMCalPtTaskTrackSelectionESD():
		AliEMCalPtTaskVTrackSelection(),
		fTrackCuts(NULL)
{
	/*
	 * Default constructor
	 */
}

//______________________________________________________________________________
AliEMCalPtTaskTrackSelectionESD::AliEMCalPtTaskTrackSelectionESD(AliESDtrackCuts* cuts):
		AliEMCalPtTaskVTrackSelection(),
		fTrackCuts(cuts)
{
	/*
	 * Constructor with cuts
	 */
}

//______________________________________________________________________________
AliEMCalPtTaskTrackSelectionESD::AliEMCalPtTaskTrackSelectionESD(
		const AliEMCalPtTaskTrackSelectionESD& ref):
		AliEMCalPtTaskVTrackSelection(ref),
		fTrackCuts(NULL)
{
	/*
	 * Copy constructor, creating a new cut object
	 */
	if(ref.fTrackCuts) fTrackCuts = new AliESDtrackCuts(*(ref.fTrackCuts));
}

//______________________________________________________________________________
AliEMCalPtTaskTrackSelectionESD& AliEMCalPtTaskTrackSelectionESD::operator=(
		const AliEMCalPtTaskTrackSelectionESD& ref)
{
	/*
	 * Assignment operator
	 */
	AliEMCalPtTaskVTrackSelection::operator=(ref);
	if(&ref != this){
		this->~AliEMCalPtTaskTrackSelectionESD();
	}
	if(ref.fTrackCuts) fTrackCuts = new AliESDtrackCuts(*(ref.fTrackCuts));
	return *this;
}

//______________________________________________________________________________
AliEMCalPtTaskTrackSelectionESD::~AliEMCalPtTaskTrackSelectionESD() {
	/*
	 * Destructor, deleting track cuts
	 */
	if(fTrackCuts) delete fTrackCuts;
}

//______________________________________________________________________________
TObjArray* AliEMCalPtTaskTrackSelectionESD::GetAcceptedTracks(
		const TClonesArray* const tracks) {
	/*
	 * Select tracks from a TClonesArray of input tracks
	 *
	 * @param tracks: TClonesArray of tracks (must not be null)
	 * @return: TObjArray of selected tracks
	 */
	if(!fListOfTracks) fListOfTracks = new TObjArray;
	else fListOfTracks->Clear();
	if(!fTrackCuts){
		AliError("Track cuts not provided");
		return fListOfTracks;
	}
	TIter trackIter(tracks);
	AliESDtrack *track(NULL);
	while((track = dynamic_cast<AliESDtrack *>(trackIter()))){
		if(fTrackCuts->AcceptTrack(track)) fListOfTracks->AddLast(track);
	}
	return fListOfTracks;
}



//______________________________________________________________________________
TObjArray* AliEMCalPtTaskTrackSelectionESD::GetAcceptedTracks(const AliVEvent* const event) {
	/*
	 * Select tracks from a virtual event
	 *
	 * @param event: AliESDEvent, via interface of virtual event (must not be null)
	 * @return: TObjArray of selected tracks
	 */
	if(!fListOfTracks) fListOfTracks = new TObjArray;
	else fListOfTracks->Clear();
	if(!fTrackCuts){
		AliError("Track cuts not provided");
		return fListOfTracks;
	}
	const AliESDEvent *esd = dynamic_cast<const AliESDEvent *>(event);
	if(!esd){
		AliError("Event not of type AliESDEvent");
		return fListOfTracks;
	}
	std::auto_ptr<TObjArray> accepted(fTrackCuts->GetAcceptedTracks(esd));
	TIter trackIter(accepted.get());
	AliESDtrack *track(NULL);
	while((track = dynamic_cast<AliESDtrack *>(trackIter()))){
		fListOfTracks->AddLast(track);
	}
	return fListOfTracks;
}

} /* namespace EMCalTriggerPtAnalysis */
