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
 * Implementation of track selection in case the analysis runs on AODs
 * For the moment it uses the AliESDtrackCuts and converts AOD tracks to
 * ESD tracks, which might change in the future when an AOD track selection
 * framework becomes available.
 *
 * Author:
 *   Markus Fasel
 */
#include <TClonesArray.h>
#include <TObjArray.h>

#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <AliESDtrack.h>
#include <AliEMCalPtTaskTrackSelectionAOD.h>

ClassImp(EMCalTriggerPtAnalysis::AliEMCalPtTaskTrackSelectionAOD)

namespace EMCalTriggerPtAnalysis {

	//______________________________________________________________________________
	AliEMCalPtTaskTrackSelectionAOD::AliEMCalPtTaskTrackSelectionAOD() :
		AliEMCalPtTaskVTrackSelection(),
		fTrackCuts(NULL),
		fFilterBits(0)
	{
		/*
		 * Main constructor
		 */
	}

	//______________________________________________________________________________
	AliEMCalPtTaskTrackSelectionAOD::AliEMCalPtTaskTrackSelectionAOD(AliESDtrackCuts* cuts, UInt_t filterbits):
		AliEMCalPtTaskVTrackSelection(),
		fTrackCuts(cuts),
		fFilterBits(filterbits)
	{
		/*
		 * Main Constructor, initalising also track cuts
		 *
		 * @param cuts: Inital track cut object
		 */
	}

	//______________________________________________________________________________
	AliEMCalPtTaskTrackSelectionAOD::AliEMCalPtTaskTrackSelectionAOD(const AliEMCalPtTaskTrackSelectionAOD& ref) :
		AliEMCalPtTaskVTrackSelection(ref),
		fTrackCuts(NULL),
		fFilterBits(ref.fFilterBits)
	{
		/*
		 * copy constructor
		 *
		 * @param ref: AOD track selection as basis for the copy
		 */
		if(ref.fTrackCuts) fTrackCuts = new AliESDtrackCuts(*(ref.fTrackCuts));
	}

	//______________________________________________________________________________
	AliEMCalPtTaskTrackSelectionAOD& AliEMCalPtTaskTrackSelectionAOD::operator=(const AliEMCalPtTaskTrackSelectionAOD& ref) {
		/*
		 * Assignment operator
		 *
		 * @param ref: AOD track selection as basis for the copy
		 * @return: reference to this cut object
		 */
		AliEMCalPtTaskVTrackSelection::operator=(ref);
		if(this != &ref){
			if(fTrackCuts) {
				delete fTrackCuts;
				fTrackCuts = NULL;
			}
			if(ref.fTrackCuts) fTrackCuts = new AliESDtrackCuts(*(ref.fTrackCuts));
		}
		return *this;
	}

	//______________________________________________________________________________
	AliEMCalPtTaskTrackSelectionAOD::~AliEMCalPtTaskTrackSelectionAOD() {
		/*
		 * Destructor, removes the track cuts
		 */
		if(fTrackCuts) delete fTrackCuts;
	}

	//______________________________________________________________________________
	TObjArray* AliEMCalPtTaskTrackSelectionAOD::GetAcceptedTracks(const TClonesArray* const tracks) {
		/*
		 * Select tracks from a list (TClonesArray) of tracks. Internally, the tracks are converted
		 * to ESD tracks and processed by the underlying AliESDtrackCut object
		 *
		 * @param tracks: TClonesArray of input tracks, under which we select the appropriate ones
		 * @return: TObjArray of selected tracks
		 */
		if(!fListOfTracks) fListOfTracks = new TObjArray;
		else fListOfTracks->Clear();
		TIter trackIter(tracks);
		AliAODTrack *track(NULL);
		while((track = dynamic_cast<AliAODTrack *>(trackIter()))){
			// First check filter bits
			if(fFilterBits && !track->TestFilterBit(fFilterBits)) continue;
			if(fTrackCuts){
				AliESDtrack copyTrack(track);
				if(fTrackCuts->AcceptTrack(&copyTrack)) fListOfTracks->AddLast(track);
			}
		}
		return fListOfTracks;
	}

	//______________________________________________________________________________
	TObjArray* AliEMCalPtTaskTrackSelectionAOD::GetAcceptedTracks(const AliVEvent* const event) {
		/*
		 * Select tracks from a list (TClonesArray) of tracks. Internally, the tracks are converted
		 * to ESD tracks and processed by the underlying AliESDtrackCut object
		 *
		 * @param tracks: TClonesArray of input tracks, under which we select the appropriate ones
		 * @return: TObjArray of selected tracks
		 */
		if(!fListOfTracks) fListOfTracks = new TObjArray;
		else fListOfTracks->Clear();
		const AliAODEvent *aod = dynamic_cast<const AliAODEvent *>(event);
		if(!aod){
			AliError("Event not of type AliAODEvent");
			return fListOfTracks;
		}
		AliAODTrack *track(NULL);
		for(int itrk = 0; itrk < event->GetNumberOfTracks(); itrk++){
			track = static_cast<AliAODTrack *>(event->GetTrack(itrk));
			// First check filter bits
			if(fFilterBits && !track->TestFilterBit(fFilterBits)) continue;
			if(fTrackCuts){
				AliESDtrack copyTrack(track);
				if(fTrackCuts->AcceptTrack(&copyTrack)) fListOfTracks->AddLast(track);
			}
		}
		return fListOfTracks;
	}

} /* namespace EMCalTriggerPtAnalysis */

