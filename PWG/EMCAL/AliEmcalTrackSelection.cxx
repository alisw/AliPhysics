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
	fListOfCuts(NULL),
	fSelectionModeAny(kFALSE)
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
AliEmcalTrackSelection& AliEmcalTrackSelection::operator=(const AliEmcalTrackSelection& ref) {
	TObject::operator=(ref);
	if(this != &ref){
		this->~AliEmcalTrackSelection();
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
AliEmcalTrackSelection::~AliEmcalTrackSelection() {
	if(fListOfTracks) delete fListOfTracks;
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
