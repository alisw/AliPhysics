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
#include <AliEmcalTrackSelection.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <AliVTrack.h>
#include <AliVEvent.h>
#include "AliVCuts.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTrackSelection)
/// \endcond

/**
 * Default consturctor, initialising objects with NULL,
 * sets acception mode to ALL
 */
AliEmcalTrackSelection::AliEmcalTrackSelection() :
	TObject(),
	fListOfTracks(NULL),
  fListOfTrackBitmaps(NULL),
  fTrackBitmap(64),
	fListOfCuts(NULL),
	fSelectionModeAny(kFALSE)
{
}

/**
 * Copy constructor, performing a flat copy
 * \param ref
 */
AliEmcalTrackSelection::AliEmcalTrackSelection(const AliEmcalTrackSelection& ref):
	TObject(ref),
	fListOfTracks(NULL),
	fListOfTrackBitmaps(NULL),
	fTrackBitmap(64),
	fListOfCuts(NULL),
	fSelectionModeAny(kFALSE)
{
	if(ref.fListOfTracks) fListOfTracks = new TObjArray(*(ref.fListOfTracks));
	if(ref.fListOfTrackBitmaps) fListOfTrackBitmaps = new TObjArray(*(ref.fListOfTrackBitmaps));
	if(ref.fListOfCuts){
	  fListOfCuts = new TObjArray;
	  fListOfCuts->SetOwner(false);
	  for(TIter cutIter = TIter(ref.fListOfCuts).Begin(); cutIter != TIter::End(); ++cutIter)
	    fListOfCuts->Add(*cutIter);
	}
}

/**
 * Assingment operator, makes a flat copy
 * \param ref Reference for the copy
 * \return Result of the copy
 */
AliEmcalTrackSelection& AliEmcalTrackSelection::operator=(const AliEmcalTrackSelection& ref) {
	TObject::operator=(ref);
	if(this != &ref){
		this->~AliEmcalTrackSelection();
		if(ref.fListOfTracks) fListOfTracks = new TObjArray(*(ref.fListOfTracks));
		if(ref.fListOfTrackBitmaps) fListOfTrackBitmaps = new TObjArray(*(ref.fListOfTrackBitmaps));
		if(ref.fListOfCuts){
		  fListOfCuts = new TObjArray;
		  fListOfCuts->SetOwner(false);
		  for(TIter cutIter = TIter(ref.fListOfCuts).Begin(); cutIter != TIter::End(); ++cutIter)
		    fListOfCuts->Add(*cutIter);
		} else fListOfCuts = NULL;
	}
	return *this;
}

/**
 * Destructor, deletes track and track cut arrays
 * In case the object has ownership over the track cuts itself, it also deletes those
 */
AliEmcalTrackSelection::~AliEmcalTrackSelection() {
	if(fListOfTracks) delete fListOfTracks;
	if(fListOfTrackBitmaps) delete fListOfTrackBitmaps;
	if(fListOfCuts) delete fListOfCuts;
}

/**
 * Add new track cuts to the list of cuts. Takes ownership over the cuts
 * \param cuts New cuts to add
 */
void AliEmcalTrackSelection::AddTrackCuts(AliVCuts *cuts){
  if(!fListOfCuts){
    fListOfCuts = new TObjArray;
    fListOfCuts->SetOwner(true);
  }
  fListOfCuts->Add(cuts);
}

/**
 * Add new track cuts to the list of cuts. Takes ownership over the cuts
 * \param cuts New cuts to add
 */
void AliEmcalTrackSelection::AddTrackCuts(TObjArray *cuts){
  TIter next(cuts);
  AliVCuts* item = 0;
  while ((item = static_cast<AliVCuts*>(next())))
  {
    AddTrackCuts(item);
  }
}

/**
 * Get the number of cut objects assigned.
 * \return The number of cut objects
 */
Int_t AliEmcalTrackSelection::GetNumberOfCutObjects() const {
  if(!fListOfCuts) return 0;
  return fListOfCuts->GetEntries();
}

/**
 * Access to track cuts at a given position
 * \param icut Cut at position in array
 * \return The cuts (NULL for invalid positions)
 */
AliVCuts* AliEmcalTrackSelection::GetTrackCuts(Int_t icut) {
  if(!fListOfCuts) return NULL;
  if(icut < fListOfCuts->GetEntries())
    return static_cast<AliVCuts *>(fListOfCuts->At(icut));
  return NULL;
}

/**
 * Select tracks from a TClonesArray of input tracks
 *
 * \param tracks TClonesArray of tracks (must not be null)
 * \return: TObjArray of selected tracks
 */
TObjArray* AliEmcalTrackSelection::GetAcceptedTracks(const TClonesArray* const tracks)
{
  if (!fListOfTracks) {
    fListOfTracks = new TObjArray;
  }
  else {
    fListOfTracks->Clear();
  }

  if (!fListOfTrackBitmaps) {
    fListOfTrackBitmaps = new TObjArray;
    fListOfTrackBitmaps->SetOwner(kTRUE);
  }
  else {
    fListOfTrackBitmaps->Clear();
  }

  TIter next(tracks);
  AliVTrack* track = 0;
  while((track = static_cast<AliVTrack*>(next()))) {
    if (IsTrackAccepted(track)) {
      fListOfTracks->AddLast(track);
    }
    else {
      fListOfTracks->AddLast(0);
    }
    fListOfTrackBitmaps->Add(new TBits(fTrackBitmap));
  }
  return fListOfTracks;
}

/**
 * Select tracks from a virtual event. Delegates selection process to function IsTrackAccepted
 *
 * \param event AliESDEvent, via interface of virtual event (must not be null)
 * \return TObjArray of selected tracks
 */
TObjArray* AliEmcalTrackSelection::GetAcceptedTracks(const AliVEvent* const event)
{
  if (!fListOfTracks) {
    fListOfTracks = new TObjArray;
  }
  else {
    fListOfTracks->Clear();
  }

  if (!fListOfTrackBitmaps) {
    fListOfTrackBitmaps = new TObjArray;
    fListOfTrackBitmaps->SetOwner(kTRUE);
  }
  else {
    fListOfTrackBitmaps->Clear();
  }

  for(int itrk = 0; itrk < event->GetNumberOfTracks(); itrk++){
    AliVTrack *trk = static_cast<AliVTrack*>(event->GetTrack(itrk));
    if (IsTrackAccepted(trk)) {
      fListOfTracks->AddLast(trk);
    }
    else {
      fListOfTracks->AddLast(trk);
    }
    fListOfTrackBitmaps->Add(new TBits(fTrackBitmap));
  }
  return fListOfTracks;
}
