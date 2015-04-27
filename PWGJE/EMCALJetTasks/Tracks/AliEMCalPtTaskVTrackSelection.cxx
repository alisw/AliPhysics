/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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
#include <TObjArray.h>
#include "AliVCuts.h"
#include "AliEMCalPtTaskVTrackSelection.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * Default consturctor, initialising objects with NULL aa
 */
AliEMCalPtTaskVTrackSelection::AliEMCalPtTaskVTrackSelection() :
	TObject(),
	fListOfTracks(NULL),
	fListOfCuts(NULL)
{
}

/**
 * Copy constructor, performing a flat copy
 * \param ref
 */
AliEMCalPtTaskVTrackSelection::AliEMCalPtTaskVTrackSelection(const AliEMCalPtTaskVTrackSelection& ref):
	TObject(ref),
	fListOfTracks(NULL),
	fListOfCuts(NULL)
{
	if(ref.fListOfTracks) fListOfTracks = new TObjArray(*(ref.fListOfTracks));
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
AliEMCalPtTaskVTrackSelection& AliEMCalPtTaskVTrackSelection::operator=(const AliEMCalPtTaskVTrackSelection& ref) {
	TObject::operator=(ref);
	if(this != &ref){
		this->~AliEMCalPtTaskVTrackSelection();
		if(ref.fListOfTracks) fListOfTracks = new TObjArray(*(ref.fListOfTracks));
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
AliEMCalPtTaskVTrackSelection::~AliEMCalPtTaskVTrackSelection() {
	if(fListOfTracks) delete fListOfTracks;
	if(fListOfCuts) delete fListOfCuts;
}

/**
 * Add new track cuts to the list of cuts. Takes ownership over the cuts
 * \param cuts New cuts to add
 */
void AliEMCalPtTaskVTrackSelection::AddTrackCuts(AliVCuts *cuts){
  if(!fListOfCuts){
    fListOfCuts = new TObjArray;
    fListOfCuts->SetOwner(true);
  }
  fListOfCuts->Add(cuts);
}

/**
 * Get the number of cut objects assigned.
 * \return The number of cut objects
 */
Int_t AliEMCalPtTaskVTrackSelection::GetNumberOfCutObjects() const {
  if(!fListOfCuts) return 0;
  return fListOfCuts->GetEntries();
}

/**
 * Access to track cuts at a given position
 * \param icut Cut at position in array
 * \return The cuts (NULL for invalid positions)
 */
AliVCuts* AliEMCalPtTaskVTrackSelection::GetTrackCuts(Int_t icut) {
  if(!fListOfCuts) return NULL;
  if(icut < fListOfCuts->GetEntries())
    return static_cast<AliVCuts *>(fListOfCuts->At(icut));
  return NULL;
}

} /* namespace EMCalTriggerPtAnalysis */
